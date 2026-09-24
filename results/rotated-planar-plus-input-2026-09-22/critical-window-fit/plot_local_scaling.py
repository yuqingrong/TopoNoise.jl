"""Render the saved near-critical fits and the unchanged full-scan counts."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import numpy as np

RUN = Path(__file__).resolve().parent
RAW = RUN.parent / "plus_input_logical_errors-raw.csv"
REPORT = json.loads((RUN / "local_fits.json").read_text())
assert hashlib.sha256(RAW.read_bytes()).hexdigest() == REPORT["raw_sha256"]
assert hashlib.sha256((RUN / "fit_local_scaling.py").read_bytes()).hexdigest() == REPORT["script_sha256"]
FITS = {tuple(key.split("/")): value for key, value in REPORT["primary_fits"].items()}
DISTANCES = (9, 11, 13, 15)
COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
PANELS = (("as", "x_only"), ("as", "z_only"),
          ("bp", "x_only"), ("bp", "z_only"))

with RAW.open(newline="") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 816
series = {}
for panel in PANELS:
    for distance in DISTANCES:
        selected = sorted((r for r in rows
                           if (r["construction"], r["error_channel"]) == panel
                           and int(r["distance"]) == distance),
                          key=lambda r: float(r["p_x"]) + float(r["p_z"]))
        assert len(selected) == 51
        assert all(r["logical_state"] == ("plus" if panel[0] == "as" else "zero")
                   and r["logical_observable"] == ("logical_x" if panel[1] == "x_only" else "logical_z")
                   and int(r["shots"]) == 50000 for r in selected)
        p = np.array([float(r["p_x"]) + float(r["p_z"]) for r in selected])
        y = np.array([int(r["logical_failures"]) / int(r["shots"]) for r in selected])
        se = np.sqrt(y * (1-y) / 50000)
        assert np.allclose(p, np.arange(51) * 0.002, rtol=0, atol=1e-14)
        assert np.allclose(y, [float(r["logical_failure_rate"]) for r in selected], rtol=0, atol=1e-14)
        assert np.allclose(se, [float(r["logical_failure_standard_error"]) for r in selected], rtol=0, atol=1e-14)
        series[panel, distance] = (p, y, se)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 10,
    "axes.edgecolor": "#404040", "axes.linewidth": 0.8,
    "grid.color": "#dddddd", "grid.linewidth": 0.6,
    "legend.fontsize": 9, "legend.framealpha": 1,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
})


def panel_title(panel):
    name = "As" if panel[0] == "as" else "Bp"
    state = "+" if panel[0] == "as" else "0"
    pauli = "X" if panel[1] == "x_only" else "Z"
    return f"{name} construction  |  " + rf"$|{state}_L\rangle$" + f"  |  {pauli}-only noise"


def fit_text(fit):
    return (rf"$p_c={fit['p_c']:.5f}\pm{fit['p_c_standard_error']:.5f}$" + "\n"
            + rf"$\nu={fit['nu']:.2f}\pm{fit['nu_standard_error']:.2f}$")


def local_data(panel, distance):
    p, y, se = series[panel, distance]
    lo, hi = FITS[panel]["window"]
    selected = (p >= lo - 1e-12) & (p <= hi + 1e-12)
    assert selected.sum() == 7
    return p[selected], y[selected], se[selected]


def collapse(ax, panel, small=False):
    fit = FITS[panel]
    all_x = []
    for distance, color in zip(DISTANCES, COLORS):
        p, y, se = local_data(panel, distance)
        x = (p - fit["p_c"]) * distance**(1/fit["nu"])
        all_x.extend(x)
        ax.errorbar(x, y, yerr=se, color=color, marker="o",
                    markersize=2 if small else 4, linestyle="none",
                    elinewidth=0.5 if small else 0.8, capsize=1 if small else 2,
                    label=f"d = {distance}", zorder=3)
    xx = np.linspace(min(all_x), max(all_x), 250)
    ax.plot(xx, np.polynomial.polynomial.polyval(xx, fit["coefficients"]),
            "--", color="#333333", linewidth=1 if small else 1.5, label="Quadratic fit")
    ax.set_xlabel(r"$(p-p_c)d^{1/\nu}$", fontsize=8 if small else 11,
                  labelpad=1 if small else 4)
    ax.set_ylabel("Logical error rate", fontsize=8 if small else 10,
                  labelpad=1 if small else 4)
    if small:
        ax.xaxis.label.set_bbox({"facecolor": "white", "edgecolor": "none", "pad": 0.6})
    ax.tick_params(labelsize=7 if small else 9, length=2 if small else 3.5)
    ax.xaxis.set_major_locator(MaxNLocator(3 if small else 5))
    ax.yaxis.set_major_locator(MaxNLocator(3 if small else 5))
    ax.grid(True, alpha=0.6)


def draw_full_panel(ax, panel):
    ax.set_title(panel_title(panel), pad=10)
    for distance, color in zip(DISTANCES, COLORS):
        p, y, se = series[panel, distance]
        ax.errorbar(p, y, yerr=se, color=color, marker="o", markersize=3,
                    linewidth=1.4, elinewidth=0.6, capsize=1.4,
                    label=f"d = {distance}", zorder=3)
    pauli = "X" if panel[1] == "x_only" else "Z"
    ax.set(xlim=(-0.002, 0.102), ylim=(-0.012, 0.53),
           xlabel="Error probability p per gate operand", ylabel=f"Logical {pauli} error rate")
    ax.set_xticks(np.arange(0, 0.101, 0.02))
    ax.set_yticks(np.arange(0, 0.51, 0.1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.grid(True, alpha=0.7)
    ax.legend(loc="lower right")
    if panel not in FITS:
        ax.text(0.17, 0.27, "No reliable threshold fit in this scan.",
                transform=ax.transAxes, fontsize=9, color="#444444")
        return
    fit = FITS[panel]
    lo, hi = fit["window"]
    ax.axvspan(lo, hi, color="#dce9e4", alpha=0.65, zorder=1)
    ax.axvline(fit["p_c"], color="#333333", linewidth=1, linestyle="--", zorder=2)
    ax.text(0.97, 0.31, f"Local fit: {lo:.3f} ≤ p ≤ {hi:.3f}\n" + fit_text(fit),
            ha="right", va="bottom", transform=ax.transAxes, fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 2})
    inset = ax.inset_axes([0.115, 0.62, 0.37, 0.30], zorder=10)
    inset.set_facecolor("white")
    collapse(inset, panel, small=True)


def save(fig, stem):
    for extension in ("png", "svg", "pdf"):
        path = RUN / f"{stem}.{extension}"
        fig.savefig(path, dpi=200)
        assert path.stat().st_size > 1000
        print(path, flush=True)
    plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize=(13, 9.4))
fig.subplots_adjust(left=0.072, right=0.98, bottom=0.105, top=0.88,
                    wspace=0.24, hspace=0.34)
fig.suptitle("Matched inputs and gates: As / Bp logical X and Z errors", y=0.975,
             fontsize=16, fontweight="bold")
fig.text(0.5, 0.93, "Near-critical fits  |  d = 9, 11, 13, 15  |  50,000 shots / point",
         ha="center", fontsize=11)
for ax, panel in zip(axes.flat, PANELS):
    draw_full_panel(ax, panel)
fig.text(0.5, 0.032,
         "Shaded bands and inset points: fitted windows only; quadratic joint scaling fits with 1,000 count-bootstrap replicates.\n"
         "Error bars and fit ± values: one statistical standard error; window/model and finite-size effects are additional.\n"
         "Ideal inputs: As all-plus, Bp all-zero. Independent noise after matched H/CNOT operands; perfect final syndrome.",
         ha="center", va="center", fontsize=8, color="#555555")
save(fig, "plus_input_logical_errors-local-fit")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.subplots_adjust(left=0.083, right=0.975, bottom=0.115, top=0.87,
                    wspace=0.24, hspace=0.38)
fig.suptitle("Near-critical scaling: As / Z and Bp / X", y=0.97,
             fontsize=16, fontweight="bold")
fig.text(0.5, 0.925, r"$P_L=A+Bx+Cx^2,\quad x=(p-p_c)d^{1/\nu}$"
         + "   |   Joint fit to all four distances", ha="center", fontsize=11)
for column, panel in enumerate((("as", "z_only"), ("bp", "x_only"))):
    fit = FITS[panel]
    lo, hi = fit["window"]
    ax = axes[0, column]
    ax.set_title(panel_title(panel), pad=10)
    for distance, color in zip(DISTANCES, COLORS):
        p, y, se = local_data(panel, distance)
        ax.errorbar(p, y, yerr=se, color=color, marker="o", markersize=4,
                    linestyle="none", elinewidth=0.8, capsize=2, label=f"d = {distance}", zorder=3)
        pp = np.linspace(lo, hi, 200)
        xx = (pp-fit["p_c"]) * distance**(1/fit["nu"])
        ax.plot(pp, np.polynomial.polynomial.polyval(xx, fit["coefficients"]),
                color=color, linewidth=1.2)
    ax.axvline(fit["p_c"], linestyle="--", color="#333333", linewidth=1)
    pauli = "Z" if panel[1] == "z_only" else "X"
    ax.set(xlim=(lo-0.0004, hi+0.0004), xlabel="Error probability p per gate operand",
           ylabel=f"Logical {pauli} error rate")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    ax.grid(True, alpha=0.7)
    ax.legend(loc="upper left", ncol=2)
    ax.text(0.97, 0.06, fit_text(fit) + "\n"
            + rf"$\chi^2/\mathrm{{dof}}={fit['reduced_chi_square']:.2f}$",
            ha="right", va="bottom", transform=ax.transAxes, fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 2})
    bottom = axes[1, column]
    bottom.set_title(f"Data collapse: {lo:.3f} ≤ p ≤ {hi:.3f}", pad=10)
    collapse(bottom, panel)
    bottom.legend(loc="upper left", fontsize=8)
fig.text(0.5, 0.032,
         "28 points and 5 fitted parameters per panel; error bars and fit ± values are one statistical standard error.\n"
         "1,000/1,000 bootstrap fits converge in each panel. Estimates are conditional on the window and model.\n"
         "Only d = 9–15 are sampled; distance and window sensitivity are documented separately. Ideal native inputs and fully matched noisy H/CNOT gates.",
         ha="center", va="center", fontsize=8.5, color="#555555")
save(fig, "near_critical_scaling")
print("Validated all 816 original counts and both 28-point fit selections; saved 6 exports.")
