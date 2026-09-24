"""Verify input provenance, independent dual curves, fits, and exports."""
import csv
import hashlib
import json
from pathlib import Path
import re
import tomllib
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image
from scipy.stats import chi2

RUN = Path(__file__).resolve().parent
REPO = RUN.parent.parent
FIT = RUN / "critical-window-fit"
RAW = RUN / "plus_input_logical_errors-raw.csv"
metadata = tomllib.loads((RUN / "metadata.toml").read_text())
assert metadata["completed_utc"] and metadata["source_unchanged_during_run"]
for relative, digest in metadata["source_sha256"].items():
    assert hashlib.sha256((RUN / "source_snapshot" / relative).read_bytes()).hexdigest() == digest
    assert hashlib.sha256((REPO / relative).read_bytes()).hexdigest() == digest

def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

rows = read_csv(RAW)
archive_dir = REPO / "results/rotated-planar-native-logical-errors-2026-09-22"
archive = read_csv(archive_dir / "native_logical_errors-raw.csv")
assert len(rows) == 816 and sum(int(r["shots"]) for r in rows) == 40_800_000
assert [r for r in rows if r["construction"] == "bp"] == [
    r for r in archive if r["construction"] == "bp"]
assert metadata["new_shots"] == metadata["reused_shots"] == 20_400_000

def key(row):
    return (row["construction"], row["error_channel"], int(row["distance"]),
            round(float(row["p_x"]) + float(row["p_z"]), 10))

indexed = {key(row): row for row in rows}
assert len(indexed) == len(rows)
duality = []
for as_channel, bp_channel in (("z_only", "x_only"), ("x_only", "z_only")):
    standardized = []
    for distance in (9, 11, 13, 15):
        for p in np.arange(51) * 0.002:
            a = indexed["as", as_channel, distance, round(p, 10)]
            b = indexed["bp", bp_channel, distance, round(p, 10)]
            n1, n2 = int(a["shots"]), int(b["shots"])
            k1, k2 = int(a["logical_failures"]), int(b["logical_failures"])
            pooled = (k1 + k2) / (n1 + n2)
            # The normal/chi-square approximation is used only where both
            # pooled expected counts per experiment are at least 5.
            if min(n1, n2) * min(pooled, 1-pooled) < 5:
                continue
            variance = pooled * (1-pooled) * (1/n1 + 1/n2)
            standardized.append((k1/n1 - k2/n2) / np.sqrt(variance))
    statistic = float(np.dot(standardized, standardized))
    dof = len(standardized)
    duality.append({"as_channel": as_channel, "bp_channel": bp_channel,
                    "points_used": dof, "chi_square": statistic,
                    "chi_square_per_point": statistic/dof,
                    "approximate_agreement_p_value": float(chi2.sf(statistic, dof)),
                    "max_absolute_standardized_difference": float(np.max(np.abs(standardized)))})

report = json.loads((FIT / "local_fits.json").read_text())
assert hashlib.sha256(RAW.read_bytes()).hexdigest() == report["raw_sha256"]
assert hashlib.sha256((FIT / "fit_local_scaling.py").read_bytes()).hexdigest() == report["script_sha256"]
bootstrap = read_csv(FIT / "bootstrap_samples.csv")
assert len(bootstrap) == 2000
for name, fit in report["primary_fits"].items():
    construction, channel = name.split("/")
    assert fit["window"] == [0.020, 0.032]
    selected = [r for r in rows if r["construction"] == construction
                and r["error_channel"] == channel
                and 0.020-1e-12 <= float(r["p_x"])+float(r["p_z"]) <= 0.032+1e-12]
    p = np.array([float(r["p_x"])+float(r["p_z"]) for r in selected])
    d = np.array([int(r["distance"]) for r in selected])
    y = np.array([int(r["logical_failures"])/int(r["shots"]) for r in selected])
    x = (p-fit["p_c"]) * d**(1/fit["nu"])
    expected = np.polynomial.polynomial.polyval(x, fit["coefficients"])
    statistic = np.sum((y-expected)**2 / (y*(1-y)/50000))
    assert len(selected) == fit["points"] == 28
    assert np.isclose(statistic, fit["chi_square"], rtol=1e-12)
    samples = [r for r in bootstrap if r["construction"] == construction and r["error_channel"] == channel]
    assert len(samples) == fit["bootstrap_successes"] == fit["bootstrap_replicates"] == 1000
    assert {int(r["replicate"]) for r in samples} == set(range(1000))
    for parameter in ("p_c", "nu"):
        values = np.array([float(r[parameter]) for r in samples])
        assert np.isclose(values.std(ddof=1), fit[parameter+"_standard_error"])
        assert np.allclose(np.quantile(values, [.025,.975]), fit[parameter+"_interval_95"])
assert len(read_csv(FIT / "window_sensitivity.csv")) == 16
assert len(read_csv(FIT / "distance_sensitivity.csv")) == 4

for stem in ("plus_input_logical_errors-local-fit", "near_critical_scaling"):
    with Image.open(FIT / (stem+".png")) as im:
        im.verify()
    ET.parse(FIT / (stem+".svg"))
    pdf = (FIT / (stem+".pdf")).read_bytes()
    assert pdf.startswith(b"%PDF-") and pdf.rstrip().endswith(b"%%EOF")
log = (RUN / "tests.log").read_text()
assert "Test Failed" not in log and "Error During Test" not in log and "ERROR:" not in log
counts = [int(m[1]) for line in log.splitlines()
          if (m := re.search(r"\|\s+(\d+)\s+(\d+)\s+(?:\d|[.])", line))]
assert len(counts) == 34 and sum(counts) == 6464
a, b = report["primary_fits"]["as/z_only"], report["primary_fits"]["bp/x_only"]
delta = a["p_c"] - b["p_c"]
se = float(np.hypot(a["p_c_standard_error"], b["p_c_standard_error"]))
output = {"raw_sha256": report["raw_sha256"], "original_counts_verified": 816,
          "independent_new_as_shots": 20_400_000, "verified_reused_bp_shots": 20_400_000,
          "test_assertions_passed": sum(counts), "testsets_passed": len(counts),
          "exports_verified": 6, "dual_curve_comparisons": duality,
          "threshold_difference_as_minus_bp": delta,
          "threshold_difference_statistical_standard_error": se,
          "threshold_difference_in_statistical_se": delta/se}
(RUN / "validation.json").write_text(json.dumps(output, indent=2)+"\n")
print(json.dumps(output, indent=2))
