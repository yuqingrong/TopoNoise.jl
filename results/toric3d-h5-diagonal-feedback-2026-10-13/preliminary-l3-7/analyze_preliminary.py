"""Apply the unchanged minimum-size guard to the user-requested diagnostics."""
import hashlib
import json
import sys
import tomllib
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parents[2]
sys.path.insert(0, str(repo / "examples"))
import analyze_toric3d_vertex_projector as analysis

selection_path = here / "selection.json"
selection = json.loads(selection_path.read_text())
raw = here / "production/toric3d-h5-vertex-projector-feedback-raw.csv"
metadata_path = here / "production/toric3d-h5-vertex-projector-feedback-metadata.toml"
meta = tomllib.loads(metadata_path.read_text())
assert selection["status"] == "diagnostic_only"
assert meta["stage"] == "preliminary_three_size" and meta["complete"]
assert meta["invalid_recoveries"] == 0 and meta["shots"] == 1000
assert meta["selection_sha256"] == analysis.sha(selection_path)
assert meta["raw_csv_sha256"] == analysis.sha(raw)
for name in ("seed", "decoder_settings", "source_hashes"):
    assert meta[name] == selection[name]
assert all(analysis.sha(repo / p) == h for p, h in meta["source_hashes"].items())
rows = analysis.read_points(raw)
expected = {(ch, *dims, p) for ch, g in selection["channels"].items()
            for dims in g["dimensions_T_Lu_Lv"] for p in g["rates"]}
assert len(rows) == len(expected) == 195
assert {analysis.key(r) for r in rows} == expected
assert {r["T"] for r in rows} == {3, 5, 7}

fits = []
for channel in analysis.CHANNELS:
    for protocol in analysis.PROTOCOLS:
        candidates = []
        for index, window in enumerate(selection["channels"][channel]["fit_windows"][protocol]):
            result = analysis.fit_window(rows, channel, protocol, window,
                                         bootstrap=1000, seed=1234)
            assert "need at least four sizes" in result["diagnostic"]
            assert result["candidate"] is None and result["pc"] is None
            assert result["bootstrap"]["accepted"] == 0
            assert result["bootstrap"]["samples"] == []
            result["candidate_index"] = index
            candidates.append(result)
        fits.append(dict(channel=channel, protocol=protocol,
                         status="diagnostic_candidates", pc=None, nu=None,
                         pc_interval=None, candidates=candidates))
report = dict(fits=fits, sensitivity=[], selection_rule=selection["rule_version"],
              scope="User-confirmed L3,5,7 preliminary diagnostics",
              user_amendment=selection["user_amendment"],
              bootstrap_execution="Not run: minimum four-size guard fails before fitting",
              sensitivity_execution="Not run: no eligible scaling fit with three sizes")
output = here / "analysis"
analysis.write_report(report, rows, output, dict(
    selection=analysis.sha(selection_path), raw_csv=analysis.sha(raw),
    metadata=analysis.sha(metadata_path), analysis=analysis.sha(analysis.__file__),
    preliminary_driver=analysis.sha(__file__),
    likelihood_engine=analysis.sha(repo / "examples/fit_toric3d_diagonal_feedback.py")))
p = output / "report.md"
p.write_text(p.read_text().replace(
    "# Direct-projector H+5 guarded analysis",
    "# Preliminary H+5 diagnostics: L=3,5,7\n\n"
    "All threshold estimates are withheld. The unchanged minimum-four-size guard "
    "rejects every candidate before fitting. The configured 1000 bootstrap "
    "replicates were therefore not executed; empty bootstrap/sensitivity files "
    "do not represent successful fits. Blank scaling panels indicate that no "
    "scaling fit was performed. Fixed-size curves and paired effects remain valid."))
print("PASS: all eight arms threshold-ineligible; no fitted thresholds or bootstrap claims")
