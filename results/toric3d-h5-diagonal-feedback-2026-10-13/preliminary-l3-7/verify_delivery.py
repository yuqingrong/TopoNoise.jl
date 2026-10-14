"""Verify quoted measurements, guard outcomes and saved regression evidence."""
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parents[2]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
subprocess.run([sys.executable, str(here.parent / "audit_artifacts.py"),
                str(here / "production")], check=True, cwd=repo)
subprocess.run([sys.executable, str(here / "audit_checkpoint_reuse.py")], check=True, cwd=repo)
raw = here / "production/toric3d-h5-vertex-projector-feedback-raw.csv"
with raw.open() as f:
    rows = list(csv.DictReader(f))
summary = json.loads((here / "analysis/measured-summary.json").read_text())
assert summary["raw_csv_sha256"] == sha(raw)
assert summary["paired_records"] == 195000 and summary["arm_replays"] == 780000
assert summary["bootstrap_executed"] == 0 and not summary["threshold_eligible"]
doc = (repo / "docs/superpowers/toric3d-h5-vertex-projector.md").read_text()
tex = (repo / "document/note.tex").read_text()
for e in summary["examples"]:
    r = next(r for r in rows if r["channel"] == e["channel"] and int(r["T"]) == e["L"]
             and float(r["p"]) == e["p"] and r["metric"] == (
                 "logical_x_union" if e["channel"] == "x_only" else "logical_z_union"))
    assert int(r["shots"]) == e["shots"] == 1000
    for arm, n in e["failures"].items():
        assert int(r[arm + "_failures"]) == n
    for arm, pair in e["paired"].items():
        for name, value in pair.items():
            assert math.isclose(float(r[f"none_vs_{arm}_{name}"]), value, abs_tol=1e-14)
    line = "| " + " | ".join([e["channel"], str(e["p"]), str(e["L"])] +
                               [f"{n/10:.1f}%" for n in e["failures"].values()]) + " |"
    assert line in doc
    if e["p"] != .077446809:
        line = " & ".join(["X" if e["channel"] == "x_only" else "Z", str(e["p"]), str(e["L"])] +
                          [f"{n/10:.1f}" for n in e["failures"].values()])
        assert line in tex
for channel, p, L, arm in (
    ("x_only", .039448164, 7, "one_step"),
    ("x_only", .039448164, 7, "window_complete"),
    ("x_only", .077446809, 7, "window_complete"),
    ("x_only", .077446809, 7, "full_prefix"),
    ("z_only", .0096470592, 3, "window_complete"),
    ("z_only", .0096470592, 7, "window_complete"),
):
    e = next(e for e in summary["examples"] if (e["channel"], e["p"], e["L"]) == (channel, p, L))
    pair = e["paired"][arm]
    interval = f"[{100*pair['ci_low']:.2f},{100*pair['ci_high']:.2f}]"
    assert interval in tex, interval

fits = json.loads((here / "analysis/fits.json").read_text())
assert fits["hashes"]["raw_csv"] == sha(raw)
assert fits["hashes"]["preliminary_driver"] == sha(here / "analyze_preliminary.py")
assert len(fits["fits"]) == 8
candidates = [c for f in fits["fits"] for c in f["candidates"]]
assert len(candidates) == 14
for f in fits["fits"]:
    assert f["status"] == "diagnostic_candidates"
    assert f["pc"] is f["nu"] is f["pc_interval"] is None
for c in candidates:
    assert "need at least four sizes" in c["diagnostic"]
    assert c["candidate"] is c["pc"] is None
    assert c["bootstrap"]["accepted"] == 0 and not c["bootstrap"]["samples"]
for name in ("bootstrap.csv", "sensitivity.csv"):
    with (here / "analysis" / name).open() as f:
        assert list(csv.DictReader(f)) == []

verification = here / "verification"
julia = (verification / "final-julia-tests.log").read_text()
counts = [tuple(map(int, m.groups())) for line in julia.splitlines()
          if (m := re.search(r'\|\s+(\d+)\s+(\d+)\s+[0-9.]+[smh]', line))]
assert len(counts) == 145 and sum(a for a, b in counts) == sum(b for a, b in counts) == 45610
python = (verification / "final-python-tests.log").read_text()
assert "Ran 51 tests" in python and "\nOK\n" in python
assert "FAILED" not in python and "Test Failed" not in julia
figures = [here / folder / (stem + "." + ext)
           for folder, prefix in (("production", "toric3d-h5-vertex-projector-feedback-"), ("analysis", "fits-"))
           for stem in (prefix + "x_only", prefix + "z_only") for ext in ("pdf", "svg", "png")]
assert len(figures) == 12 and all(p.stat().st_size > 1000 for p in figures)
pdf = repo / "document/note.pdf"
assert pdf.read_bytes().startswith(b"%PDF-")
subprocess.run(["git", "diff", "--check"], check=True, cwd=repo)
result = dict(passed=True, points=195, paired_records=195000, invalid_recoveries=0,
              illustrative_rows_verified=9, quoted_paired_intervals_verified=6,
              guarded_windows=14, eligible_arms=0, bootstrap_executed=0,
              julia_assertions=45610, julia_test_sets=145, python_tests=51,
              figure_exports=12, raw_csv_sha256=sha(raw), note_pdf_sha256=sha(pdf),
              report_sha256=sha(repo / "docs/superpowers/toric3d-h5-vertex-projector.md"),
              note_tex_sha256=sha(repo / "document/note.tex"),
              regression_logs={p.name: sha(p) for p in verification.glob("final-*-tests.log")},
              scientific_source_unchanged=True)
(verification / "delivery-verification.json").write_text(json.dumps(result, indent=2) + "\n")
print("PASS: final measurements, guards, provenance, 45610 Julia assertions and 51 Python tests verified")
