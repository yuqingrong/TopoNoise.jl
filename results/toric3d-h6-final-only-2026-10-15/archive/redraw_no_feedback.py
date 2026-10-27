"""Redraw the saved no-feedback fits in the paper theme; no fitting or sampling."""
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/toponoise-h6-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
COLORS = {3: '#7F7F7F', 5: '#0072B2', 7: '#E69F00', 9: '#009E73',
          11: '#CC79A7', 13: '#0072B2', 17: '#E69F00', 21: '#CC79A7', 25: '#7F7F7F'}


def style_axis(axis, inset=False):
    axis.grid(False, which='both')
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color('0.25')
        spine.set_linewidth(.8)
    if inset:
        axis.tick_params(labelsize=8, length=3, pad=2)
        axis.xaxis.label.set_fontsize(9)
        axis.yaxis.label.set_fontsize(9)
        axis.xaxis.labelpad = 1
        axis.yaxis.labelpad = 2


def draw_points(axis, rows, fit=None):
    for size in sorted({row['steps'] for row in rows}):
        selected = sorted((row for row in rows if row['steps'] == size), key=lambda row: row['p'])
        x = np.array([row['p'] for row in selected])
        if fit is not None:
            x = (x - fit['pc']) * size ** (1 / fit['nu'])
        y = np.array([row['rate'] for row in selected])
        errors = np.maximum(0, [y - np.array([row['ci_low'] for row in selected]),
                               np.array([row['ci_high'] for row in selected]) - y])
        axis.errorbar(x, y, yerr=errors, marker='o', color=COLORS[size], label=fr'$L={size}$',
                      ms=2.6 if fit else 3.8, lw=.9 if fit else 1.5,
                      elinewidth=.65 if fit else .85, capsize=1.4 if fit else 2.2, zorder=3)


def main():
    report = json.loads((STUDY / 'data/threshold-report.json').read_text())
    with (STUDY / 'data/threshold-points.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key in ('steps', 'shots', 'failures'):
            row[key] = int(row[key])
        for key in ('p', 'rate', 'ci_low', 'ci_high'):
            row[key] = float(row[key])
        assert row['rate'] == row['failures'] / row['shots']
        assert row['ci_low'] <= row['rate'] <= row['ci_high']

    assert report['fits']['x_only']['status'] == 'diagnostic_only'
    assert report['fits']['x_only']['pc'] is None
    assert report['fits']['z_only']['status'] == 'accepted'
    with plt.style.context(HERE / 'paper-no-grid.mplstyle'):
        figure, axes = plt.subplots(1, 2, figsize=(13.2, 6.2), sharey=True)
        for axis, channel in zip(axes, ('x_only', 'z_only')):
            result = report['fits'][channel]
            fit = result['candidate']
            selected = sorted((r for r in rows if r['channel'] == channel), key=lambda r: (r['steps'], r['p']))
            # Reconstruct saved predictions, proving the plot uses the original fitted model.
            scaled = np.array([(r['p'] - fit['pc']) / fit['p_scale'] *
                               (r['steps'] / fit['L_reference']) ** (1 / fit['nu']) for r in selected])
            predictions = expit(np.polynomial.polynomial.polyval(scaled, fit['coefficients']))
            np.testing.assert_allclose(predictions, fit['fitted_probabilities'], rtol=1e-12, atol=1e-14)
            raw_file = STUDY / 'data' / ('production-x' if channel == 'x_only' else 'production-z') / 'raw.csv'
            with raw_file.open() as stream:
                raw = {(int(r['steps']), float(r['p'])): (int(r['shots']), int(r['failures']))
                       for r in csv.DictReader(stream)}
            assert raw == {(r['steps'], r['p']): (r['shots'], r['failures']) for r in selected}
            draw_points(axis, selected)
            style_axis(axis)
            axis.set(title='X faults · BP+OSD' if channel == 'x_only' else 'Z faults · MWPM',
                     xlabel=fr'Physical {channel[0].upper()} fault probability $p$', ylim=(0, .46))
            rates = np.array(result['rates'])
            span = rates[-1] - rates[0]
            axis.set_xlim(rates[0] - .05 * span, rates[-1] + .05 * span)
            axis.set_xticks([.002, .004, .006, .008] if channel == 'x_only'
                            else [.0045, .0050, .0055, .0060, .0065])
            axis.tick_params(axis='x', labelsize=10)
            axis.legend(loc='upper right', ncol=1, fontsize=9)
            inset = axis.inset_axes([.12, .66, .48, .26])
            draw_points(inset, selected, fit=fit)
            x = (np.array([r['p'] for r in selected]) - fit['pc']) * np.array([r['steps'] for r in selected]) ** (1 / fit['nu'])
            grid = np.linspace(x.min(), x.max(), 300)
            normalized = grid / (fit['p_scale'] * fit['L_reference'] ** (1 / fit['nu']))
            inset.plot(grid, expit(np.polynomial.polynomial.polyval(normalized, fit['coefficients'])),
                       color='0.2', lw=1, zorder=2)
            inset.set(xlabel=r'$(p-p_0)L^{1/\nu_0}$' if channel == 'x_only' else r'$(p-p_c)L^{1/\nu}$',
                      ylabel='Logical error')
            inset.set_title('Diagnostic collapse · rejected' if channel == 'x_only' else 'Accepted scaling collapse',
                            fontsize=8, pad=4)
            style_axis(inset, inset=True)
            inset.locator_params(axis='x', nbins=4)
            inset.locator_params(axis='y', nbins=3)
            inset.margins(x=.05, y=.12)
            if channel == 'x_only':
                annotation = ('Diagnostic fit only\n' + fr'$p_0={fit["pc"]:.6f}$' + '\n' +
                              fr'$\nu_0={fit["nu"]:.3f}$' + '\nNo accepted threshold\n' +
                              r'Fit $p$-value: $4.50\times10^{-7}$')
            else:
                low, high = result['bootstrap']['pc_interval']
                axis.axvspan(low, high, ymax=.58, color='0.2', alpha=.08, lw=0, zorder=0)
                axis.axvline(result['pc'], ymax=.58, color='0.35', ls='--', lw=1, zorder=1)
                annotation = ('Accepted threshold fit\n' + fr'$p_c={fit["pc"]:.6f}$ ({100*fit["pc"]:.4f}%)' + '\n' +
                              fr'$\nu={fit["nu"]:.3f}$' + '\n' +
                              f'95% CI: {100*low:.4f}%–{100*high:.4f}%')
            axis.text(.04, .49, annotation, transform=axis.transAxes, va='top', fontsize=10,
                      bbox=dict(facecolor='white', edgecolor='none', alpha=.96, pad=2), zorder=5)
        axes[0].set_ylabel('Logical state failure rate')
        figure.suptitle('H+6 · No intermediate feedback · One final recovery', y=.975, fontsize=14)
        figure.text(.075, .062, r'T = L · 95% Wilson error bars · X: 1,000 shots/point, ideal $|00\rangle$ seed · '
                    r'Z: 10,000 shots/point, noisy-H $|++\rangle$ seed', fontsize=9)
        figure.text(.075, .025, 'Z shaded band: conditional 95% bootstrap interval. X fitted candidate fails the threshold checks.', fontsize=9)
        figure.subplots_adjust(left=.075, right=.98, bottom=.19, top=.85, wspace=.20)
        figure.canvas.draw()
        for axis in [*axes, *(child for ax in axes for child in ax.child_axes)]:
            assert not any(line.get_visible() for line in axis.get_xgridlines() + axis.get_ygridlines())
            assert all(spine.get_visible() for spine in axis.spines.values())
        output = STUDY / 'figures/h6-threshold-failure-rates.pdf'
        figure.savefig(output)
        plt.close(figure)
    print('Saved PDF in paper theme; 80 observed points and both saved models verified; all grids off.')


if __name__ == '__main__':
    main()
