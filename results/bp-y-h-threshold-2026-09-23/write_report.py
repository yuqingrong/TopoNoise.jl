"""Write the results report from validated counts and fit summaries."""
import json
from pathlib import Path

from analysis import CHANNELS, load_points

RUN = Path(__file__).resolve().parent


def main():
    validation = json.loads((RUN / "validation.json").read_text())
    if validation["status"] != "passed":
        raise ValueError("Report requires a successful artifact validation")
    fits = json.loads((RUN / "local_fits.json").read_text())["primary_fits"]
    selection = json.loads((RUN / "selection.json").read_text())
    points = load_points()
    names = {"y_only": "Y-only", "hadamard": "Exact H=(X+Z)/√2"}
    table = ["| Noise | Logical component | p_c | ν | χ²/dof | Fit window |",
             "|---|---|---|---|---|---|"]
    intervals, sensitivity_lines = [], []
    sensitivity = json.loads((RUN / "sensitivity.json").read_text())
    for channel in CHANNELS:
        for component in ("x", "z"):
            key = f"{channel}/{component}"
            f = fits[key]
            if f["status"] == "success":
                table.append(f"| {names[channel]} | {component.upper()} | "
                             f"{f['p_c']:.6f} ± {f['p_c_standard_error']:.6f} | "
                             f"{f['nu']:.3f} ± {f['nu_standard_error']:.3f} | "
                             f"{f['reduced_chi_square']:.2f} | {f['window'][0]:.3f}–{f['window'][1]:.3f} |")
                low, high = f["p_c_interval_95"]
                nu_low, nu_high = f["nu_interval_95"]
                intervals.append(f"- {names[channel]}/{component.upper()}: 95% statistical interval for p_c "
                                 f"[{low:.6f}, {high:.6f}] and for ν [{nu_low:.3f}, {nu_high:.3f}]; "
                                 f"{f['bootstrap_successes']}/{f['bootstrap_replicates']} "
                                 f"bootstrap fits succeeded; goodness-of-fit p = {f['goodness_of_fit_p']:.3f}.")
                variants = [s for s in sensitivity if s["panel"] == key and "p_c" in s]
                quadratic = [s for s in variants if s["kind"] == "window" and s["degree"] == 2]
                cubic = [s for s in variants if s["kind"] == "window" and s["degree"] == 3]
                distance = [s for s in variants if s["kind"] == "distance"]
                for label, group in (("quadratic windows", quadratic), ("cubic windows", cubic)):
                    if group:
                        sensitivity_lines.append(f"- {names[channel]}/{component.upper()}, {label}: "
                                                 f"p_c={min(s['p_c'] for s in group):.6f}–"
                                                 f"{max(s['p_c'] for s in group):.6f}; "
                                                 f"χ²/dof={min(s['reduced_chi_square'] for s in group):.2f}–"
                                                 f"{max(s['reduced_chi_square'] for s in group):.2f}.")
                for s in distance:
                    sensitivity_lines.append(f"- {names[channel]}/{component.upper()}, distances "
                                             f"{', '.join(map(str,s['distances']))}: p_c={s['p_c']:.6f}, "
                                             f"ν={s['nu']:.3f}, χ²/dof={s['reduced_chi_square']:.2f}.")
            else:
                status = "No supported crossing in the scanned range" if f["status"] == "unresolved" else "Local fit unreliable"
                table.append(f"| {names[channel]} | {component.upper()} | {status} | — | — | — |")
    pilot_points = [p for p in points if p["stage"] == "pilot"]
    production_points = [p for p in points if p["stage"] == "production"]
    refined = sum(p["shots"] == 10000 for p in pilot_points)
    lines = [
        "# Bp circuit under Y-only and exact Hadamard faults", "",
        "![Logical error curves and local threshold fits](logical_errors.png)", "",
        "[Figure PDF](logical_errors.pdf) · [Near-critical fits and collapse](near_critical_scaling.pdf) · "
        "[Raw counts](raw_counts.csv) · [Fit summary](fit_summary.csv)", "",
        "## Results", "", *table, "",
        "The ± values are one statistical standard error. Each logical-X and logical-Z rate is a "
        "component marginal: a joint logical-Y outcome contributes to both. The logical-Z curves "
        "do not support a positive finite-size crossing in the scanned range. This observation "
        "does not by itself establish a zero asymptotic threshold.", "", *intervals, "",
        "These are finite-size encoding-channel estimates conditional on this Bp schedule, "
        "the perfect final syndrome round, the independent CSS matching decoder, and the selected "
        "scaling window. Statistical intervals exclude window, polynomial, and distance-selection effects.", "",
        "## Exact physical model and logical readout", "",
        "Each active data-qubit operand receives an independent post-gate fault with total probability p. "
        "The two channels are E_Y(ρ)=(1−p)ρ+pYρY and E_H(ρ)=(1−p)ρ+pHρH, "
        "where H=(X+Z)/√2. Each Y event flips both Pauli-frame components together. "
        "H is replayed as an actual Hadamard in a full Stim stabilizer trajectory, retaining its "
        "quantum action through all subsequent gates and measurements.", "",
        "The data use the current native Bp gate schedule, boundary orientation x_ns, at distances "
        "9, 11, 13, and 15. Fault locations are the stored H gates and individual CNOT operands. "
        "Input preparation, the reference qubit, and final "
        "measurements are ideal; there is no idle or measurement noise.", "",
        "To characterize both logical components with one fixed noisy schedule, the ideal input is "
        "(U_Bp†⊗I_R)|Φ_L+⟩. Final checks are followed by the commuting observables Z_L Z_R and X_L X_R. "
        "Their signs are XORed with the X- and Z-sector corrections predicted by the existing "
        "uniform-weight PyMatching decoder, with correlations disabled. Only final syndromes enter "
        "the decoder. The two failure bits f_X,f_Z give joint counts n00,n10,n01,n11, with "
        "P_X=(n10+n11)/N and P_Z=(n01+n11)/N. This Bell-reference experiment characterizes the "
        "logical channel; the reference=0 branch reproduces the native all-zero Bp input.", "",
        "## Sampling and fitting", "",
        f"The archive contains {validation['points']} points and {validation['shots']:,} shots: "
        f"{len(pilot_points)} pilot points ({sum(p['shots'] for p in pilot_points):,} shots) and "
        f"{len(production_points)} fresh production points ({sum(p['shots'] for p in production_points):,} shots). "
        f"The pilot uses p=0:0.005:0.1 and 2,000 shots per point; {refined} ambiguous-region points "
        "were extended to 10,000 shots. Production uses 50,000 independent shots per point and Δp=0.002. "
        "Both logical marginals are retained from each trial.", "",
        "Crossing candidates exclude p=0 and pilot logical rates outside (0.005,0.45). "
        "The selected adjacent-distance crossings must form a compact cluster, with supported "
        "d=9 versus d=15 ordering reversal. A broad but supported extreme-distance reversal "
        "triggers further pilot sampling instead of an immediate no-crossing conclusion. "
        "The full numerical rules and candidate lists are in analysis.py and selection.json. "
        "The initial Y selection was retained unchanged while the H pilot was refined; no H "
        "production data existed before its window was frozen. Both selection records are archived.", "",
        "Production centers are " + ", ".join(f"{names[c]}/{component.upper()}: {choice['center']:.3f}"
             for c in CHANNELS for component in ("x", "z")
             if (choice := selection["panels"][f"{c}/{component}"])["status"] == "selected") + ". "
        "Each production window extends ±0.010; the primary fit uses ±0.006.", "",
        "The local model is P_L=A+B*x+C*x², x=(p−p_c)d^(1/ν), with inverse observed-binomial-variance "
        "weights and four optimizer starts. ν is searched over (0.3,5). "
        "The [PanQEC threshold-fitting tutorial](https://panqec.readthedocs.io/en/latest/tutorials/Computing%20threshold.html) "
        "describes this near-transition expansion. All four distances are fitted jointly. "
        "Each of 1,000 multinomial bootstrap draws resamples joint outcome counts and supplies "
        "both logical marginals, preserving their correlation. Variances are recomputed each draw. "
        "Accepted fits require interior, identifiable parameters, physical increasing probabilities, "
        "at least 95% bootstrap success, and goodness-of-fit p≥0.01.", "",
        "## Sensitivity", "", *sensitivity_lines, "",
        "All attempted quadratic/cubic fits at halfwidths 0.004,0.006,0.008,0.010 and both distance "
        "omissions are retained in sensitivity.csv and sensitivity.json, including rejected or "
        "poorly fitting cases. These alternatives share samples and are not independent estimates.", "",
        "## Verification and provenance", "",
        "The standalone tests verify noiseless Bell checks at d=3,9,11,13,15; 73 independent "
        "Yao statevector histories at d=3; Y histories against the original Julia Pauli propagation "
        "and decoder; exact H batch replay; multiple faults and H²=I; probability endpoints; "
        "seeded repetition; checkpoint continuation; changed-source rejection; a synthetic known "
        "threshold/exponent; row-order invariance; saturation/flat-fit rejection; and correlated "
        "bootstrap behavior. The dispersed H pilot is retained as a regression fixture for refinement.", "",
        "The complete Julia package regression suite passed on retry after granting access to its "
        "normal plotting cache; the original sandbox denial and successful retry logs are both retained. "
        "Exact test exit codes and log hashes are in test_status.json. validation.json verifies "
        "every checkpoint count, CSV marginal and error bar, source hash, production window, "
        "fit objective, bootstrap interval, plotted series, and all six figure exports.", "",
        "Both final PNG figures were visually inspected for readable labels, complete panels, "
        "error bars, fit-window shading, scaling insets, and caption spacing. The PDF and SVG "
        "exports use the same figure objects. Styling follows the September 22 native logical-error "
        "figures: the same distance colors, boxed axes, per-panel legends, upper-left collapse "
        "insets, and typography. figure_style_revision.json records the reference and unchanged "
        "data/fit hashes for this presentation-only revision.", "",
        "metadata.json records runtime versions, seeds, circuit/sampler fingerprints, and completion. "
        "source_snapshot retains the actual working-tree Julia sources and dependency manifests. "
        "Sampling uses deterministic seeds per 1,000-shot chunk and atomic checkpoints, independent "
        "of worker scheduling. Reproducibility is conditional on the recorded Stim version and machine "
        "architecture. No existing package API or earlier result archive was modified.", "",
        "## Reproduce", "",
        "Run from the repository root in the installed project environment. Completed points resume "
        "without resampling; a changed sampler or environment fingerprint is rejected. "
        "requirements-python.txt records the exact Python package versions used.", "",
        "```sh",
        "MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\",
        "  .CondaPkg/.pixi/envs/default/bin/python \\",
        "  results/bp-y-h-threshold-2026-09-23/run_experiment.py --stage all --workers 4",
        "",
        "MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \\",
        "  results/bp-y-h-threshold-2026-09-23/analysis.py",
        "MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \\",
        "  results/bp-y-h-threshold-2026-09-23/plot_results.py",
        "MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \\",
        "  results/bp-y-h-threshold-2026-09-23/validate_results.py",
        "```", "",
        "The sampler reads archived circuits.toml directly. export_circuits.jl regenerates schedules "
        "and Yao fixtures from matching Julia sources and refuses to overwrite a different export.", "",
        "For a fresh replay, make a new directory under results/ and copy the Python scripts, "
        "circuits.toml, source_hashes.toml, source_snapshot/, oracle_fixtures.toml, fixtures/, and "
        "requirements-python.txt into it. Leave checkpoints/, metadata.json, and selection*.json "
        "uncopied. Run the new directory's runner, analysis, and plot scripts with the same commands. "
        "The unchanged root seed, chunk size, architecture, and recorded package versions reproduce "
        "the simulated counts. Final validation additionally requires that directory's completed "
        "test logs and test_status.json, as in this archive.", "",
    ]
    (RUN / "README.md").write_text("\n".join(lines))
    print("Wrote README.md from validated results.")


if __name__ == "__main__":
    main()
