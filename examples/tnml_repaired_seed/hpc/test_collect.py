"""Synthetic reporting checks; no Julia or scheduler jobs are launched."""
import copy
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


C = module("hpc_collect_test_target", "collect.py")


def toml(document):
    flat = ''.join(f'{json.dumps(k)} = {json.dumps(v)}\n' for k, v in document.items() if not isinstance(v, dict))
    return flat + ''.join(f'\n[{k}]\n' + toml(v) for k, v in document.items() if isinstance(v, dict))


def fixture(root, manifest, task, *, count=None, invalid=None, first_failure=False, raw_variant=False):
    path = root / task["output"]
    (path / "shots").mkdir(parents=True, exist_ok=True)
    (path / "exact").mkdir(exist_ok=True)
    N = task["shots"]
    count = N if count is None else count
    run_id = f"synthetic-{task['stage']}-{task['task_id']}-chi{task['chi']}"
    hashes = manifest["scientific_source_sha256"]
    metadata = dict(L=task["L"], p=task["p"], shots=N, chi=task["chi"], cutoff=task["cutoff"],
                    rng_seed=task["seed"], exact=task["exact"], negative_tolerance=manifest["negative_tolerance"],
                    tn_uses_both_records=True, run_id=run_id,
                    project_sha256=hashes["Project.toml"],
                    source_sha256={k:v for k,v in hashes.items() if k not in ("Project.toml", "Manifest.toml")})
    metadata.update(julia_version="1.12.4", ITensors_version="0.9.10", ITensorMPS_version="0.3.24",
                    PythonCall_version="0.9.31", python_version="3.12.9", python_executable="/shared/python",
                    pymatching_version="2.3.1", julia_threads=1, blas_threads=1,
                    TopoNoise_version="0.1.0", TopoNoise_project_path="/shared/project",
                    active_project_path="/shared/project/Project.toml", physical_noise="synthetic Z faults",
                    preparation="synthetic seed repair", tn_conditioning="seed plus final",
                    mwpm_conditioning="final", rng_algorithm="Random.Xoshiro")
    if "Manifest.toml" in hashes:
        metadata["manifest_sha256"] = hashes["Manifest.toml"]
    (path / "metadata.toml").write_text(toml(metadata))
    records, failures, valid = [], 0, 0
    for i in range(1, N+1):
        truth = int(first_failure and i == 1)
        record = dict(shot_index=i, raw_fault_count=8, raw_fault_indices=[1] if truth or raw_variant else [],
                      seed_syndrome=[False, False], seed_recovery=[False, False], final_syndrome=[False, False],
                      recovery_mask=0, raw_seed_mask=truth, transfer_mask=0, true_mask=truth, mwpm_failure=False)
        record["record_id"] = C.record_identity(record)
        records.append(record)
        if i > count:
            continue
        okay = i != invalid
        shot = dict(L=task["L"], p=task["p"], chi=task["chi"], cutoff=task["cutoff"], shot_index=i,
                    run_id=run_id, record_id=record["record_id"], true_mask=truth, mwpm_failure=False,
                    status="complete" if okay else "invalid", valid=okay, logical_mask=0,
                    TNfailure=bool(truth), probabilities=[1.,0.,0.,0.], elapsed_seconds=1.)
        if task["exact"] and okay:
            exact = dict(shot)
            exact.pop("TNfailure")
            (path / "exact" / f"shot_{i:04d}.toml").write_text(toml(exact))
            shot["exact_comparison"] = dict(ml_tie_aware_match=True, ml_strict_match=True, posterior_linf=0.)
        (path / "shots" / f"shot_{i:04d}.toml").write_text(toml(shot))
        valid += okay
        failures += truth if okay else 0
    (path / "input-corpus.toml").write_text(toml(dict(run_id=run_id, requested_shots=N, rng_seed=task["seed"]))
        + ''.join('\n[[records]]\n' + toml(r) for r in records))
    (path / "summary.toml").write_text(toml(dict(L=task["L"], p=task["p"], run_id=run_id,
        requested_shots=N, valid_shots=valid, unknown_shots=N-valid, tn_failures=failures)))
    (path / "supervision.json").write_text(json.dumps(dict(status="completed", exit_code=0)))
    return path


class CollectIntegrity(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="tnml-hpc-collect-synthetic-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest = C.HPC.create_manifest(self.root)
        self.production = C.HPC.stage_tasks(self.manifest, "production")

    def test_absent_whole_shards_stay_in_planned_denominator(self):
        task = self.production[0]
        fixture(self.root, self.manifest, task, first_failure=True)
        report = C.collect_report(self.root)
        point = next(p for p in report["production_points"] if p["p"] == task["p"])
        self.assertEqual((point["requested_shots"], point["validated_shots"], point["unknown_shots"]), (1000,50,950))
        self.assertEqual(point["tn_failures"], 1)
        self.assertEqual(point["outcome_bounds"], [0.001, .951])
        self.assertIsNone(point["failure_rate"])
        self.assertIsNone(point["wilson95"])
        self.assertEqual(report["production_totals"]["requested_shots"], 22000)

    def test_complete_point_has_rate_only_after_all_twenty_shards(self):
        p = self.production[0]["p"]
        tasks = [task for task in self.production if task["p"] == p]
        for index, task in enumerate(tasks):
            fixture(self.root, self.manifest, task, first_failure=index == 0)
        point = next(pnt for pnt in C.collect_report(self.root)["production_points"] if pnt["p"] == p)
        self.assertTrue(point["complete"])
        self.assertEqual(point["failure_rate"], .001)
        self.assertEqual(point["wilson95"], list(C.BASE.wilson(1,1000)))

    def test_wrong_metadata_makes_entire_shard_unknown(self):
        task = self.production[0]
        path = fixture(self.root,self.manifest,task)
        metadata = path/"metadata.toml"
        metadata.write_text(metadata.read_text().replace(f'"rng_seed" = {task["seed"]}', '"rng_seed" = 1'))
        row = C.inspect_task(self.root,task,self.manifest)
        self.assertFalse(row["identity_valid"])
        self.assertEqual((row["validated_shots"],row["unknown_shots"]),(0,50))

    def test_cross_cutoff_and_source_hash_are_not_accepted(self):
        for kind in ("cutoff","source"):
            with self.subTest(kind=kind):
                task=self.production[0]
                path=fixture(self.root,self.manifest,task)
                metadata=path/"metadata.toml"
                text=metadata.read_text()
                if kind=="cutoff":
                    text=text.replace(f'"cutoff" = {task["cutoff"]}', '"cutoff" = 0.1')
                else:
                    text=text.replace(self.manifest["scientific_source_sha256"]["examples/tnml_repaired_seed/run.jl"],"wrong")
                metadata.write_text(text)
                row=C.inspect_task(self.root,task,self.manifest)
                self.assertEqual(row["validated_shots"],0)
                self.assertFalse(row["identity_valid"])

    def test_invalid_missing_exact_and_missing_output_are_unknown(self):
        task=self.production[0]
        path=fixture(self.root,self.manifest,task,count=49,invalid=3)
        (path/"exact/shot_0002.toml").unlink()
        row=C.inspect_task(self.root,task,self.manifest)
        self.assertEqual((row["validated_shots"],row["unknown_shots"]),(47,3))
        self.assertFalse(row["complete"])

    def test_exact_decision_is_recomputed_from_actual_reference(self):
        task=self.production[0]
        path=fixture(self.root,self.manifest,task)
        exact=path/"exact/shot_0001.toml"
        exact.write_text(exact.read_text().replace('"logical_mask" = 0','"logical_mask" = 1')
                         .replace('[1.0, 0.0, 0.0, 0.0]','[0.0, 1.0, 0.0, 0.0]'))
        row=C.inspect_task(self.root,task,self.manifest)
        self.assertEqual((row["validated_shots"],row["unknown_shots"]),(49,1))
        self.assertTrue(any("exact" in error["reason"] for error in row["shot_errors"]))

    def test_record_tampering_invalidates_the_corpus(self):
        task=self.production[0]
        path=fixture(self.root,self.manifest,task)
        corpus=path/"input-corpus.toml"
        corpus.write_text(corpus.read_text().replace('"raw_fault_indices" = []','"raw_fault_indices" = [1]',1))
        row=C.inspect_task(self.root,task,self.manifest)
        self.assertEqual(row["validated_shots"],0)
        self.assertTrue(any("record" in error for error in row["errors"]))

    def test_duplicate_planned_seed_and_output_are_rejected(self):
        for field in ("seed","output","task_id","shard_id"):
            with self.subTest(field=field):
                manifest=copy.deepcopy(self.manifest)
                tasks=C.HPC.stage_tasks(manifest,"production")
                tasks[1][field]=tasks[0][field]
                with self.assertRaisesRegex(ValueError,"duplicate"):
                    C.validate_plan(manifest)

    def test_paired_feasibility_requires_identical_records_not_just_seed(self):
        tasks=C.HPC.stage_tasks(self.manifest,"feasibility")
        pair_id=tasks[0]["pair_id"]
        pair=[t for t in tasks if t["pair_id"]==pair_id]
        fixture(self.root,self.manifest,pair[0])
        fixture(self.root,self.manifest,pair[1],raw_variant=True)
        report=C.collect_report(self.root)
        row=next(x for x in report["feasibility_pairs"] if x["pair_id"]==pair_id)
        self.assertFalse(row["same_raw_histories"])
        self.assertIsNone(row["posterior_linf"])

    def test_zero_and_smoke_are_separate_from_production(self):
        for phase in ("zero","smoke"):
            task=C.HPC.stage_tasks(self.manifest,phase)[0]
            fixture(self.root,self.manifest,task)
        report=C.collect_report(self.root)
        self.assertEqual(report["production_totals"]["validated_shots"],0)
        self.assertEqual(report["production_totals"]["requested_shots"],22000)
        self.assertTrue(all(point["requested_shots"]==0 and point["wilson95"] is None for point in report["analytical_zero"]))
        self.assertEqual(report["diagnostics"]["zero"]["validated_shots"],1)
        self.assertEqual(report["diagnostics"]["smoke"]["validated_shots"],2)

    def test_missing_or_failed_supervisor_never_qualifies_complete_shard(self):
        task=self.production[0]
        for supervision in (None, dict(status="completed",exit_code=2), dict(status="decode_timeout",exit_code=-15)):
            with self.subTest(supervision=supervision):
                path=fixture(self.root,self.manifest,task)
                if supervision is None:
                    (path/"supervision.json").unlink()
                else:
                    (path/"supervision.json").write_text(json.dumps(supervision))
                row=C.inspect_task(self.root,task,self.manifest)
                self.assertFalse(row["complete"])
                self.assertNotEqual(row["status"],"complete")

    def test_timeout_preserves_prefix_and_unknown_requested_histories(self):
        task=self.production[0]
        path=fixture(self.root,self.manifest,task,count=2,first_failure=True)
        (path/"supervision.json").write_text(json.dumps(dict(status="decode_timeout",exit_code=-15)))
        row=C.inspect_task(self.root,task,self.manifest)
        self.assertEqual((row["validated_shots"],row["unknown_shots"],row["tn_failures"]),(2,48,1))
        self.assertEqual(row["status"],"blocked")
        point=C.collect_report(self.root)["production_points"][0]
        self.assertEqual(point["unknown_shots"],998)
        self.assertIsNone(point["failure_rate"])

    def test_mixed_production_runtime_or_model_is_rejected(self):
        for field in ("python_version","python_executable","julia_version","physical_noise","tn_conditioning"):
            with self.subTest(field=field):
                fixture(self.root,self.manifest,self.production[0])
                path=fixture(self.root,self.manifest,self.production[1])
                metadata=C.read_toml(path/"metadata.toml")
                metadata[field]="different"
                (path/"metadata.toml").write_text(toml(metadata))
                with self.assertRaisesRegex(ValueError,"mixed production"):
                    C.collect_report(self.root)

    def test_missing_runtime_and_wrong_thread_count_make_shard_unknown(self):
        for field,value in (("python_version",None),("julia_threads",2)):
            with self.subTest(field=field):
                task=self.production[0]
                path=fixture(self.root,self.manifest,task)
                metadata=C.read_toml(path/"metadata.toml")
                if value is None:
                    del metadata[field]
                else:
                    metadata[field]=value
                (path/"metadata.toml").write_text(toml(metadata))
                row=C.inspect_task(self.root,task,self.manifest)
                self.assertFalse(row["identity_valid"])
                self.assertEqual(row["unknown_shots"],50)

    def test_matching_histories_across_chi_are_compared_despite_different_run_ids(self):
        tasks=C.HPC.stage_tasks(self.manifest,"feasibility")
        pair=[task for task in tasks if task["pair_id"]==tasks[0]["pair_id"]]
        for task in pair:
            fixture(self.root,self.manifest,task)
        report=C.collect_report(self.root)
        result=report["feasibility_pairs"][0]
        self.assertEqual(result["status"],"compared")
        self.assertTrue(result["same_raw_histories"])
        self.assertEqual(result["posterior_linf"],0.0)
        self.assertTrue(result["logical_decisions_match"])

    def test_cli_writes_only_fully_audited_production_points_and_optional_plot(self):
        p=self.production[0]["p"]
        for index,task in enumerate(task for task in self.production if task["p"]==p):
            fixture(self.root,self.manifest,task,first_failure=index==0)
        with redirect_stdout(io.StringIO()):
            result=C.main(["--root",str(self.root),"--plot"])
        saved=json.loads((self.root/"collected/report.json").read_text())
        self.assertEqual(sum(point["failure_rate"] is not None for point in saved["production_points"]),1)
        self.assertEqual(len(result["figures"]),2)
        self.assertTrue(all(Path(path).stat().st_size>1000 for path in result["figures"]))


if __name__ == "__main__":
    unittest.main()
