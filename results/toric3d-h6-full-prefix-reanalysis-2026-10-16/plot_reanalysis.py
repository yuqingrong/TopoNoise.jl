"""Draw PDF diagnostics from saved follow-up results; never fit or select data."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'examples'))
from plotting_theme import paper_theme, style_axes

MODELS = {
    'quadratic': ('Leading quadratic', '#0072B2'),
    'cubic': ('Leading cubic', '#7F7F7F'),
    'corr_w0p5': (r'Correction $L^{-1/2}$', '#E69F00'),
    'corr_w1': (r'Correction $L^{-1}$', '#009E73'),
    'corr_w2': (r'Correction $L^{-2}$', '#CC79A7'),
    'corr_linear_w1': (r'Offset + slope, $L^{-1}$', '#D55E00'),
}


def save(fig, name):
    path = HERE/'figures'/name
    path.parent.mkdir(exist_ok=True)
    for ax in fig.axes:
        style_axes(ax)
        assert not any(g.get_visible() for g in ax.get_xgridlines()+ax.get_ygridlines())
    fig.savefig(path, metadata={'Title': path.stem,
                               'Subject': 'Exploratory reanalysis of existing H+6 full-prefix samples.'})
    plt.close(fig)


def x_crossings(x):
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True, sharey=True)
    for ax, pair in zip(axes.flat, [(3, 5), (5, 7), (7, 9), (9, 11)]):
        rows = [r for r in x['crossing']['bands'] if (r['L_small'], r['L_large']) == pair]
        rows.sort(key=lambda r: r['p'])
        assert len(rows) == 7
        ps = [r['p'] for r in rows]
        ax.fill_between(ps, [r['ci_low'] for r in rows], [r['ci_high'] for r in rows],
                        color='#0072B2', alpha=.16, label='Simultaneous 95% band')
        ax.plot(ps, [r['difference'] for r in rows], 'o-', color='#0072B2', label='Measured difference')
        ax.axhline(0, color='.35', lw=.9, ls='--')
        ax.set_title(fr'$L={pair[1]}$ minus $L={pair[0]}$')
        ax.set_xlabel(r'Gate-output fault probability $p$')
        ax.set_ylabel(r'Difference in failure probability')
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.margins(x=.03)
    axes[0, 0].legend(loc='upper left', fontsize=8.5)
    regions = x['crossing']['secondary_all_pairs']['components']
    text = '; '.join(f'[{a:.5f}, {b:.5f}]' for a, b in regions) or 'empty'
    fig.suptitle('X: adjacent-size differences with simultaneous uncertainty', fontsize=14)
    fig.subplots_adjust(left=.10, right=.98, bottom=.17, top=.90, hspace=.35, wspace=.25)
    fig.text(.10, .065, f'All ten size pairs: common-crossing-compatible region {text}', fontsize=9)
    fig.text(.10, .028, 'Exact binomial cell intervals with Bonferroni coverage; linear interpolation between measured probabilities.', fontsize=8.5)
    save(fig, 'x-crossing-uncertainty.pdf')


def z_stability(z):
    assert len(z['fits']) == 36
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.6), sharey=True)
    handles = []
    for ax, window in zip(axes, ['full', 'narrow']):
        for model, (label, color) in MODELS.items():
            rows = sorted((r for r in z['fits'] if r['window'] == window and r['model'] == model),
                          key=lambda r: r['min_size'])
            assert len(rows) == 3
            line, = ax.plot([r['min_size'] for r in rows], [r['pc'] for r in rows],
                           color=color, marker='o', lw=1.2, label=label)
            if window == 'full': handles.append(line)
            for r in rows:
                if r['deviance_p_value'] < .01:
                    ax.plot(r['min_size'], r['pc'], marker='o', ms=6.5, mfc='white',
                            mec=color, mew=1.4, ls='none')
                if r['bound_hits'] or not r['optimizer_success'] or not r['monotone_increasing']:
                    ax.plot(r['min_size'], r['pc'], marker='x', ms=8, mew=1.7, color='black', ls='none')
        ax.set(xlabel=r'Minimum included size $L_{\min}$', xticks=[3, 5, 7], xlim=(2.7, 7.3),
               title='Full window: 0.025 to 0.040' if window == 'full' else 'Narrow window: 0.0275 to 0.0375')
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        ax.axvspan(6.65, 7.3, color='.92', zorder=-1)
    axes[0].set_ylabel(r'Fitted critical probability $p_c$')
    fig.suptitle('Z: dependence on size cutoff, fit window and correction model', fontsize=14)
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(.52, .065), fontsize=9)
    fig.subplots_adjust(left=.09, right=.98, bottom=.28, top=.86, wspace=.15)
    fig.text(.09, .022, 'Model estimates, not confidence intervals. Open circles: poor fit. Shaded Lmin = 7: three sizes only. x: numerical check failed.', fontsize=8.5)
    save(fig, 'z-finite-size-sensitivity.pdf')


def z_profiles(z):
    assert len(z['profiles']) == 4
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True, sharey=True)
    for ax, profile in zip(axes.flat, z['profiles']):
        rows = profile['evaluations']
        assert profile['minimum_agreement']
        values = [max(0., r['lr']) if r['optimizer_success'] and r['lr'] is not None else np.nan for r in rows]
        color = '#0072B2' if profile['model'] == 'quadratic' else '#009E73'
        ax.plot([r['pc'] for r in rows], values, color=color, lw=1.5)
        ax.axhline(profile['lr_cutoff'], color='.35', ls='--', lw=.9)
        ax.plot(profile['free_pc'], 0, marker='o', color=color)
        for interval in profile['regions']:
            ax.axvspan(interval['lower'], interval['upper'], color=color, alpha=.10, zorder=-1)
        label = 'Leading quadratic' if profile['model'] == 'quadratic' else r'Quadratic + $L^{-1}$ correction'
        fit = next(f for f in z['fits'] if f['fit_id'] == profile['fit_id'])
        if fit['deviance_p_value'] < .01:
            label += ' (poor fit)'
        ax.set(title=fr'$L_{{\min}}={profile["min_size"]}$: '+label,
               xlabel=r'Fixed critical probability $p_c$', ylabel=r'$2\Delta$ negative log-likelihood',
               ylim=(0, 12), xlim=(.025, .04))
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=1))
    fig.suptitle('Z: conditional likelihood profiles', fontsize=14)
    fig.subplots_adjust(left=.09, right=.98, bottom=.14, top=.89, hspace=.35, wspace=.23)
    fig.text(.09, .035, 'Shading: approximate 95% profile region within each model. Other parameters reoptimized; model uncertainty excluded.', fontsize=8.5)
    save(fig, 'z-threshold-profiles.pdf')


if __name__ == '__main__':
    x = json.loads((HERE/'x-report.json').read_text())
    z = json.loads((HERE/'z-report.json').read_text())
    with paper_theme():
        x_crossings(x)
        z_stability(z)
        z_profiles(z)
    print('Saved three no-grid diagnostic PDFs.')
