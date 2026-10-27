#!/usr/bin/env python3
"""Plot the already-computed, frozen drop-smallest-size sensitivity; never refit."""
import argparse
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/h6-full-prefix-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'examples'))
from plotting_theme import paper_theme, save_figure, style_axes

COLORS = {3: '#7F7F7F', 5: '#0072B2', 7: '#E69F00', 9: '#009E73', 11: '#CC79A7'}


def plot(report, output):
    for channel, fit in report['fits'].items():
        sensitivity = next(v for v in fit['sensitivity'] if v['name'] == 'drop_smallest')
        candidate = sensitivity['candidate']
        show_fit = candidate is not None and not candidate['bound_hits']
        excluded = min(fit['sizes'])
        rows = [r for r in report['marginals'] if r['channel'] == channel
                and r['arm'] == 'full_prefix' and r['L'] != excluded and r['p'] in fit['rates']]
        with paper_theme():
            fig, ax = plt.subplots(figsize=(8.2, 5.5))
            inset = ax.inset_axes([.12, .69, .38, .23]) if show_fit else None
            xs = []
            for size in sorted({r['L'] for r in rows}):
                group = [r for r in rows if r['L'] == size]
                p = np.array([r['p'] for r in group])
                y = np.array([r['rate'] for r in group])
                error = np.maximum(0, [y - np.array([r['ci_low'] for r in group]),
                                      np.array([r['ci_high'] for r in group]) - y])
                ax.errorbar(p, y, yerr=error, color=COLORS[size], marker='o',
                            label=fr'$L={size}$', elinewidth=.85)
                if inset is not None:
                    x = (p - candidate['pc']) * size**(1 / candidate['nu'])
                    xs.extend(x)
                    inset.errorbar(x, y, yerr=error, color=COLORS[size], marker='o',
                                   ms=2.3, lw=.8, elinewidth=.5, capsize=1)
            if inset is not None:
                grid = np.linspace(min(xs), max(xs), 300)
                normalized = grid / (candidate['p_scale'] * candidate['L_reference']**(1 / candidate['nu']))
                inset.plot(grid, expit(np.polynomial.polynomial.polyval(normalized, candidate['coefficients'])),
                           color='0.2', ls='--', lw=1)
                style_axes(inset, inset=True)
                inset.set(xlabel=r'$(p-p_c)L^{1/\nu}$', ylabel='Logical error')
                inset.xaxis.set_major_locator(MaxNLocator(4))
                inset.set_title('Sensitivity fit: diagnostic only', fontsize=8)
            style_axes(ax)
            ax.set(title=f'H+6 · Full-prefix {channel[0].upper()} recovery · Excluding L = {excluded}',
                   xlabel=r'Gate-output fault probability $p$',
                   ylabel=f'Logical {channel[0].upper()} state failure rate', ylim=(0, 1))
            ax.margins(x=.04)
            ax.legend(loc='upper right', ncol=2, fontsize=9)
            if show_fit:
                ax.axvline(candidate['pc'], color='0.45', ls='--', lw=.85, zorder=0)
                annotation = (fr'Candidate $p_c={candidate["pc"]:.5f}$; $\nu={candidate["nu"]:.2f}$'
                              + '\n' + f'Goodness-of-fit p = {candidate["deviance_p_value"]:.3g}')
            else:
                annotation = 'Sensitivity fit unavailable' if candidate is None else 'Sensitivity fit hits parameter bounds'
            annotation += '\nDiagnostic only; no independent threshold claim'
            ax.text(.97, .13, annotation, transform=ax.transAxes, ha='right', fontsize=9)
            fig.subplots_adjust(left=.12, right=.97, bottom=.16, top=.90)
            fig.text(.12, .025, 'Frozen drop-smallest sensitivity · 95% Wilson bars · No refitting or new data selection', fontsize=8)
            save_figure(fig, output / f'{channel[0]}-full-prefix-drop-L{excluded}-fit')
            plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plot(json.loads(args.report.read_text()), args.output)
