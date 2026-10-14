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
    "do not represent successful fits. Figures contain measured curves only, "
    "without scaling panels. Fixed-size curves and paired effects remain valid.").replace(
    "Statistical bootstrap intervals are conditional on the frozen window and scaling model. "
    "Model, size, and window sensitivity is listed separately in sensitivity.csv. "
    "A diagnostic candidate is not a threshold.",
    "No bootstrap intervals or sensitivity fits were computed. The candidate entries "
    "below are frozen scout windows rejected by the size guard, not fitted thresholds."))

# Presentation-only replacement for the generic report's empty scaling panels.
# The scientific fitting engine and its minimum-four-size guard are unchanged.
import matplotlib.pyplot as plt
for channel in analysis.CHANNELS:
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharey=True,
                             constrained_layout=True)
    for ax, protocol in zip(axes, analysis.PROTOCOLS):
        for L, marker in zip((3, 5, 7), ("o", "s", "^")):
            points = sorted((r for r in rows if r["channel"] == channel and r["T"] == L),
                            key=lambda r: r["p"])
            y = [r[f"{protocol}_failures"] / r["shots"] for r in points]
            ci = [analysis.engine.wilson_interval(r[f"{protocol}_failures"], r["shots"])
                  for r in points]
            ax.errorbar([r["p"] for r in points], y,
                        yerr=[[v - c[0] for v, c in zip(y, ci)],
                              [c[1] - v for v, c in zip(y, ci)]],
                        marker=marker, markersize=3, linewidth=1, capsize=2, label=f"L={L}")
        ax.set(title=protocol.replace("_", " "), xlabel="post-gate p", ylim=(-.02, 1.02))
        ax.grid(alpha=.2)
    label = ("X-only: logical-X union / |00_L> failure" if channel == "x_only"
             else "Z-only: state-independent logical-Z union")
    axes[0].set_ylabel("logical failure probability")
    axes[0].legend()
    fig.suptitle("Preliminary H+5 — " + label + "\n"
                 "Measured curves with Wilson 95% intervals; no scaling fit (3 sizes, at least 4 required)",
                 fontsize=12)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(output / f"fits-{channel}.{ext}", dpi=180)
    plt.close(fig)
# CSV and SVG exporters use CRLF or trailing spaces; normalize text exports.
for path in [*output.glob("*.csv"), *output.glob("*.svg")]:
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()).rstrip() + "\n")
print("PASS: all eight arms threshold-ineligible; no fitted thresholds or bootstrap claims")
