"""Focused integrity tests; all records here are synthetic and never production."""
import copy
import csv
import hashlib
import io
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location("matched_analysis", Path(__file__).with_name("analyze.py"))
AN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AN)


def row(shots=1000):
    result = dict(channel="x_only", steps=5, Ly=5, Lz=5, p=.01,
                  shots=shots, target_shots=1000, seed=123, point_seed="0123456789abcdef",
                  complete=shots == 1000, elapsed_seconds=1.0)
    bins = [shots - 10, 1, 2, 3, 1, 1, 1, 1]
    for mask, count in enumerate(bins):
        result[f"joint_{mask:02}"] = count
    for index, arm in enumerate(AN.ARMS):
        k = sum(n for mask, n in enumerate(bins) if mask & (1 << index))
        result[f"{arm}_failures"] = k
        result[f"{arm}_logical_1_failures"] = k
        result[f"{arm}_logical_2_failures"] = 0
    return result


def selection():
    return dict(shots=1000, arms=list(AN.ARMS), bootstrap_replicates=200,
                bootstrap_seed=123, channels={ch: dict(sizes=[5], rates=[.01],
                fit_rates={arm: [.01] for arm in AN.ARMS}) for ch in AN.CHANNELS},
                shards={"x": dict(channel="x_only", size=5, rate=.01, seed=123),
                        "z": dict(channel="z_only", size=5, rate=.01, seed=124)})


def toml_fixture(data):
    """Small test-only TOML writer; fixtures never reach production paths."""
    lines = []
    def scalar(value):
        return json.dumps(value, separators=(",", ":"))
    def emit(table, prefix=""):
        for key,value in sorted(table.items()):
            if not isinstance(value,dict) and not (isinstance(value,list) and value and isinstance(value[0],dict)):
                lines.append(f'{json.dumps(key)} = {scalar(value)}')
        for key,value in sorted(table.items()):
            name=(prefix+"." if prefix else "")+json.dumps(key)
            if isinstance(value,dict):
                lines.extend(["",f"[{name}]"]); emit(value,name)
            elif isinstance(value,list) and value and isinstance(value[0],dict):
                for item in value:
                    lines.extend(["",f"[[{name}]]"]); emit(item,name)
    emit(data)
    return ("\n".join(lines)+"\n").encode()


def write_mirror_fixture(study):
    frozen=selection()
    runtime=dict(julia_threads=1,blas_threads=1,python_executable="/usr/bin/python3",
                 packages={key:"test" for key in ("numpy","scipy","ldpc","pymatching")},
                 environment={key:"1" for key in ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS")})
    runtime["environment"].update(JULIA_CONDAPKG_BACKEND="Null",JULIA_PYTHONCALL_EXE="/usr/bin/python3")
    frozen.update(batch_size=25,expected_runtime=runtime,frozen_at_unix=0,protocols={ch:{} for ch in AN.CHANNELS})
    r=row()
    config=dict(channel="x_only",sizes=[5],rates=[.01],shots=1000,seed=123,batch_size=25,dimensions=[[5,5,5]])
    state=dict(completed_shots=1000,elapsed_seconds=1.0,
               logical_1_counts=[r[arm+"_logical_1_failures"] for arm in AN.ARMS],
               logical_2_counts=[r[arm+"_logical_2_failures"] for arm in AN.ARMS],
               failures=[r[arm+"_failures"] for arm in AN.ARMS],joint_counts=[r[key] for key in AN.JOINT_COLUMNS])
    shared=dict(format_version=1,config=config,runtime=runtime,source_hashes={},protocol={},
                point_seeds_hex=[r["point_seed"]],generation=1,started_at_unix=1.0,updated_at_unix=2.0,
                decoder_build_seconds_by_size={"5":0.0},decoder_build_count_by_size={"5":1})
    payload=toml_fixture(dict(shared,states=[state]))
    checkpoint=f'payload_sha256 = "{hashlib.sha256(payload).hexdigest()}"\n'.encode()+payload
    stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=AN.COLUMNS,extrasaction="ignore")
    writer.writeheader();writer.writerow(r);raw=stream.getvalue().encode()
    metadata=toml_fixture(dict(shared,complete=True,expected_points=1,completed_points=1,completed_shots=1000,
                              checkpoint_sha256=AN.digest(checkpoint),raw_sha256=AN.digest(raw)))
    folder=study/"archive/runs/x";folder.mkdir(parents=True)
    for name,contents in (("checkpoint.toml",checkpoint),("metadata.toml",metadata),("raw.csv",raw)):
        (folder/name).write_bytes(contents)
    return frozen,folder


class IntegrityTests(unittest.TestCase):
    def test_counts_and_joint_marginals(self):
        self.assertEqual(AN.validate_row(row(), selection())['shots'], 1000)
        for field in ("shots", "final_only_failures", "joint_00", "one_step_logical_1_failures"):
            bad = row()
            bad[field] += 1
            with self.subTest(field=field), self.assertRaises(ValueError):
                AN.validate_row(bad, selection())

    def test_duplicate_points_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            AN.validate_grid([row(), row()], selection())

    def test_wrong_geometry_and_unselected_rate_rejected(self):
        for changes in ({"Ly": 7}, {"p": .011}, {"target_shots": 9999}, {"seed": 99}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                AN.validate_row(dict(row(), **changes), selection())

    def test_frozen_point_seed_rejected(self):
        frozen = selection()
        frozen["shards"]["x"]["point_seed_hex"] = "fedcba9876543210"
        with self.assertRaisesRegex(ValueError, "RNG seed"):
            AN.validate_row(row(), frozen)

    def test_paired_interval_orientation_and_boundaries(self):
        for positive, negative, n in [(0, 0, 100), (100, 0, 100), (0, 100, 100), (12, 5, 100)]:
            value = AN.paired_interval(positive, negative, n)
            reverse = AN.paired_interval(negative, positive, n)
            self.assertAlmostEqual(value["difference"], (positive-negative)/n)
            self.assertLessEqual(value["ci_low"], value["difference"])
            self.assertGreaterEqual(value["ci_high"], value["difference"])
            self.assertGreaterEqual(value["ci_low"], -1)
            self.assertLessEqual(value["ci_high"], 1)
            self.assertAlmostEqual(value["ci_low"], -reverse["ci_high"])
        self.assertLess(AN.paired_interval(0, 0, 100)["ci_low"], 0)
        with self.assertRaises(ValueError):
            AN.paired_interval(70, 50, 100)

    def test_joint_contrasts_use_pairing(self):
        r = row()
        contrasts = AN.paired_contrasts(r)
        self.assertEqual(len(contrasts), 3)
        for contrast in contrasts:
            a, b = contrast["reference"], contrast["comparison"]
            expected = (r[b+"_failures"]-r[a+"_failures"])/r["shots"]
            self.assertAlmostEqual(contrast["difference"], expected)

    def test_no_fit_until_entire_grid_complete(self):
        with patch.object(AN.H6, "analyze_channel", side_effect=AssertionError("premature fit")):
            self.assertEqual(AN.fit_all([row()], selection(), complete=False), {})

    def test_fit_adapter_supplies_both_real_channels_and_frozen_windows(self):
        rows = [row(), dict(row(), channel="z_only", seed=124)]
        calls = []
        def fitted(adapted, channel, **kwargs):
            self.assertEqual({r["channel"] for r in adapted}, set(AN.CHANNELS))
            self.assertEqual(kwargs["bootstrap"], 200)
            for r in adapted:
                self.assertEqual(r["failures"], rows[AN.CHANNELS.index(r["channel"])][AN.ARMS[len(calls)//2]+"_failures"])
            calls.append(channel)
            return {"status": "diagnostic_only", "pc": None, "nu": None}
        with patch.object(AN.H6, "analyze_channel", side_effect=fitted):
            result = AN.fit_all(rows, selection(), complete=True)
        self.assertEqual(len(result), 6)
        self.assertEqual(calls, list(AN.CHANNELS)*3)

    def test_payload_checksum_detects_edit(self):
        import hashlib
        payload = b'format_version = 1\ngeneration = 2\n\n[config]\nshots = 10000\n'
        marker = f'payload_sha256 = "{hashlib.sha256(payload).hexdigest()}"\n'.encode()
        content = payload.replace(b'\n\n[config]', b'\n'+marker+b'\n[config]')
        AN.validate_payload(content)
        with self.assertRaisesRegex(ValueError, "payload"):
            AN.validate_payload(content.replace(b'shots = 10000', b'shots = 10001'))

    def test_checkpoint_metadata_csv_roundtrip_and_runtime_rejection(self):
        with tempfile.TemporaryDirectory(prefix="h6-matched-test-",dir="/private/tmp") as temporary:
            study=Path(temporary);frozen,folder=write_mirror_fixture(study)
            rows,provenance,runtime=AN.load_inputs(study,frozen,{})
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]["shots"],1000)
            self.assertEqual(provenance["x"]["generation"],1)
            bad=copy.deepcopy(frozen);bad["expected_runtime"]["julia_threads"]=2
            with self.assertRaisesRegex(ValueError,"frozen runtime"):
                AN.load_inputs(study,bad,{})

    def test_active_mirror_race_is_deferred_but_stale_corruption_rejected(self):
        with tempfile.TemporaryDirectory(prefix="h6-matched-test-",dir="/private/tmp") as temporary:
            frozen,folder=write_mirror_fixture(Path(temporary))
            (folder/"raw.csv").write_text("in-progress mirror\n")
            with self.assertRaises(AN.ActiveSnapshot):
                AN.stable_mirrors(folder,retries=1)
            for file in folder.iterdir():
                os.utime(file,(1,1))
            with self.assertRaisesRegex(ValueError,"invalid checkpoint/mirrors"):
                AN.stable_mirrors(folder,retries=1)


if __name__ == "__main__":
    unittest.main()
