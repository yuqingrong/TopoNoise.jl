"""Publication exports from frozen raw counts and verified local-fit parameters."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import numpy as np

from analysis import CHANNELS, DISTANCES, model
from sampler import atomic_json

RUN = Path(__file__).resolve().parent
# Match the September 22 native/matched-lattice figures, including distance colors.
COLORS = {9: "#0072B2", 11: "#E69F00", 13: "#009E73", 15: "#CC79A7"}
NAMES = {"y_only": "Y-only noise", "hadamard": "Exact H noise"}


def selected_rows(rows, channel, component, distance, stage=None, window=None):
    chosen = [r for r in rows if r["channel"] == channel and int(r["distance"]) == distance
              and (stage is None or r["stage"] == stage)
              and (window is None or window[0]-1e-12 <= float(r["p"]) <= window[1]+1e-12)]
    chosen.sort(key=lambda r: float(r["p"]))
    p = np.array([float(r["p"]) for r in chosen])
    y = np.array([float(r[f"logical_{component}_rate"]) for r in chosen])
    se = np.array([float(r[f"logical_{component}_standard_error"]) for r in chosen])
    return p, y, se


def save_figure(figure, name):
    for extension in ("png", "pdf", "svg"):
        figure.savefig(RUN / f"{name}.{extension}", dpi=200, facecolor="white")


def collapse(axis, rows, channel, component, fit, small=False):
    all_x = []
    for d in DISTANCES:
        p, y, se = selected_rows(rows, channel, component, d, "production", fit["window"])
        x = (p-fit["p_c"]) * d**(1/fit["nu"])
        all_x.extend(x)
        axis.errorbar(x, y, yerr=se, fmt="o-", ms=1.8 if small else 3,
                      capsize=1.0 if small else 1.4, lw=0.9 if small else 1.4,
                      elinewidth=0.6, color=COLORS[d])
    grid = np.linspace(min(all_x), max(all_x), 250)
    axis.plot(grid, np.polynomial.polynomial.polyval(grid, fit["coefficients"]),
              color="#333333", lw=1.1)
    axis.set_xlabel(r"$(p-p_c)d^{1/\nu}$", fontsize=8 if small else 10, labelpad=1,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4})
    axis.set_ylabel("logical error" if small else f"Logical {component.upper()} error rate",
                    fontsize=8 if small else 10, labelpad=1)
    axis.tick_params(labelsize=7 if small else 10, length=2 if small else 3.5)
    axis.grid(True, alpha=0.5 if small else 0.7)
    if small:
        axis.xaxis.set_major_locator(MaxNLocator(4))
        axis.yaxis.set_major_locator(MaxNLocator(3))
        axis.set_facecolor("#ffffff")
        axis.patch.set_alpha(1.0)


def main():
    fits_document = json.loads((RUN / "local_fits.json").read_text())
    raw = RUN / "raw_counts.csv"
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == fits_document["raw_sha256"]
    assert hashlib.sha256((RUN / "analysis.py").read_bytes()).hexdigest() == fits_document["script_sha256"]
    with raw.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    fits = fits_document["primary_fits"]
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 12, "axes.labelsize": 10,
        "axes.edgecolor": "#404040", "axes.linewidth": 0.8,
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.titleweight": "normal", "axes.labelcolor": "black", "text.color": "black",
        "grid.color": "#dddddd", "grid.linewidth": 0.6,
        "legend.fontsize": 9, "legend.framealpha": 1,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })
    figure, axes = plt.subplots(2, 2, figsize=(13, 9.4))
    figure.subplots_adjust(left=0.072, right=0.98, top=0.88, bottom=0.105, wspace=0.24, hspace=0.32)
    figure.suptitle(r"Bp circuit: Y-only and $H=(X+Z)/\sqrt{2}$ logical errors",
                    fontsize=16, fontweight="bold", y=0.975)
    figure.text(0.5, 0.93,
                "d = 9, 11, 13, 15   |   2,000–10,000 pilot shots   |   50,000 production shots / point",
                ha="center", fontsize=11)
    audit = []
    for i, channel in enumerate(CHANNELS):
        for j, component in enumerate(("x", "z")):
            axis = axes[i, j]
            key = f"{channel}/{component}"
            fit = fits[key]
            for d in DISTANCES:
                for stage in ("pilot", "production"):
                    p, y, se = selected_rows(rows, channel, component, d, stage)
                    if not len(p):
                        continue
                    pilot = stage == "pilot"
                    axis.errorbar(p, y, yerr=se, color=COLORS[d],
                                  fmt="o-" if pilot else "o", markersize=3,
                                  linewidth=1.4, elinewidth=0.6, capsize=1.4,
                                  label=f"d = {d}" if pilot else "_nolegend_",
                                  zorder=3 if pilot else 4)
                    audit.append(dict(panel=key, distance=d, stage=stage,
                                      p=p.tolist(), rate=y.tolist(), standard_error=se.tolist()))
            axis.set_title(f"Bp construction  |  {NAMES[channel]}  |  logical {component.upper()}", pad=10)
            axis.set_xlim(-0.002, max(0.102, max(float(r["p"]) for r in rows)+0.002))
            axis.set_ylim(-0.012, 0.53)
            axis.set_xticks(np.arange(0, 0.101, 0.02))
            axis.set_yticks(np.arange(0, 0.501, 0.1))
            axis.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            axis.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            axis.grid(True, alpha=0.7, zorder=0)
            axis.set_xlabel("Error probability p per gate operand")
            axis.set_ylabel(f"Logical {component.upper()} error rate")
            axis.legend(loc="lower right")
            if fit["status"] == "success":
                axis.axvspan(*fit["window"], color="#526271", alpha=0.08, zorder=0)
                axis.axvspan(*fit["p_c_interval_95"], color="#303b49", alpha=0.15, zorder=1)
                axis.axvline(fit["p_c"], color="#303b49", lw=1.0, ls="--", zorder=1)
                grid = np.linspace(*fit["window"], 150)
                for d in DISTANCES:
                    axis.plot(grid, model(np.array(fit["theta_internal"]), grid, d),
                              color=COLORS[d], lw=1.6, zorder=3)
                label = ("Local finite-size fit\n"
                         + rf"$p_c={fit['p_c']:.5f}\pm{fit['p_c_standard_error']:.5f}$" + "\n"
                         + rf"$\nu={fit['nu']:.2f}\pm{fit['nu_standard_error']:.2f}$")
                # The H curve remains low across the scan; keep the label above it.
                label_y = 0.69 if channel == "hadamard" else 0.31
                axis.text(0.97, label_y, label, transform=axis.transAxes,
                          ha="right", va="bottom", fontsize=9)
                inset = axis.inset_axes([0.12, 0.61, 0.37, 0.31], zorder=10)
                collapse(inset, rows, channel, component, fit, small=True)
            else:
                text = ("No supported finite-size crossing\nin the scanned p range"
                        if fit["status"] == "unresolved" else "Local scaling fit not reliable")
                axis.text(0.17, 0.27, text, transform=axis.transAxes, fontsize=9,
                          color="#444444", ha="left", va="bottom")
    figure.text(0.5, 0.026,
                "Independent Y/H faults after each gate operand; ideal Bell reference and final syndrome; CSS decoding. Seed: 20260923.\n"
                "Error bars: one binomial standard error. Light/dark shading: primary fit window / 95% threshold interval. Insets: fitted points only.\n"
                "Thresholds are finite-size estimates for this circuit and decoder; quoted errors are statistical. Joint logical Y counts in both components.",
                ha="center", fontsize=8, va="center", color="#555555")
    save_figure(figure, "logical_errors")
    plt.close(figure)

    successful = [(f"{channel}/{component}", fits[f"{channel}/{component}"])
                  for channel in CHANNELS for component in ("x", "z")
                  if fits[f"{channel}/{component}"]["status"] == "success"]
    if successful:
        zoom, axes = plt.subplots(len(successful), 2, figsize=(12.2, 3.7*len(successful)), squeeze=False)
        zoom.subplots_adjust(left=0.085, right=0.97, top=0.90, bottom=0.12, hspace=0.38, wspace=0.25)
        zoom.suptitle("Local finite-size fits and scaling collapse", fontsize=16, fontweight="bold", y=0.98)
        for i, (key, fit) in enumerate(successful):
            channel, component = key.split("/")
            left, right = axes[i]
            for d in DISTANCES:
                p, y, se = selected_rows(rows, channel, component, d, "production", fit["window"])
                left.errorbar(p, y, yerr=se, fmt="o", color=COLORS[d], ms=3, capsize=1.4,
                              elinewidth=0.6, label=f"d = {d}")
                grid = np.linspace(*fit["window"], 150)
                left.plot(grid, model(np.array(fit["theta_internal"]), grid, d), color=COLORS[d], lw=1.4)
            left.axvline(fit["p_c"], color="#526271", ls="--", lw=0.9)
            left.set_xlabel("Error probability p per gate operand")
            left.set_ylabel(f"Logical {component.upper()} error rate")
            left.set_title(f"{NAMES[channel]}: p$_c$ = {fit['p_c']:.5f} ± {fit['p_c_standard_error']:.5f}", fontsize=11)
            left.grid(True, alpha=0.7)
            left.legend(loc="upper left")
            collapse(right, rows, channel, component, fit)
            right.set_title(f"ν = {fit['nu']:.2f} ± {fit['nu_standard_error']:.2f}; "
                            f"χ²/dof = {fit['reduced_chi_square']:.2f}", fontsize=11)
        zoom.text(0.085, 0.025, "Quadratic local model; 1,000 joint-count bootstrap replicates. Statistical uncertainties only.",
                  color="#555555", fontsize=8)
    else:
        zoom, axis = plt.subplots(figsize=(9, 3))
        axis.axis("off")
        axis.text(0.5, 0.5, "No reliable local threshold fit was identified.", ha="center", va="center")
    save_figure(zoom, "near_critical_scaling")
    plt.close(zoom)
    atomic_json(RUN / "plot_data.json", dict(raw_sha256=fits_document["raw_sha256"],
                fit_sha256=hashlib.sha256((RUN / "local_fits.json").read_bytes()).hexdigest(), series=audit))
    print("Saved both figures as PNG, PDF, SVG and recorded every plotted point.")


if __name__ == "__main__":
    main()
