#!/usr/bin/env python3
"""Plot the completed targeted study; never simulate, refit, or pool data.

Defaults read analysis/report.json and write exactly five PDFs in figures/.
The two data figures use fresh full-prefix state-failure counts (either of the
two logical parities fails) and pointwise 95% Wilson intervals. Diagnostic
figures display saved analysis values and flags.
"""
from pathlib import Path
import argparse
import json
import math
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'examples'))
from plotting_theme import PALETTE, paper_theme, style_axes

MODELS = {
    'quadratic': ('Leading quadratic', '#0072B2', 'o', '-'),
    'corr_w1': (r'Quadratic + $L^{-1}$', '#009E73', 's', '--'),
}
COLORS = dict(zip((5, 7, 9, 11, 13, 15, 17), (*PALETTE, '#56B4E9')))
MARKERS = dict(zip((5, 7, 9, 11, 13, 15, 17), ('o', 's', '^', 'D', 'v', 'P', 'X')))
OUTPUTS = ('x-data.pdf', 'z-data.pdf', 'x-crossing.pdf',
           'z-size-sensitivity.pdf', 'z-profile.pdf')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def validate_report(report, plan):
    """Check completeness against the frozen job grid before creating figures."""
    require(report.get('schema_version') == 1, 'Unsupported report schema')
    expected = {}
    for job in plan['jobs']:
        for size in job['sizes']:
            for p in job['rates']:
                key = (job['channel'], size, p)
                require(key not in expected, 'Duplicate frozen job cell')
                expected[key] = job['shots']
    found = {}
    for short, channel in [('x', 'x_only'), ('z', 'z_only')]:
        data = report['fresh'][short]
        require(data['channel'] == channel and data['dataset'] == 'fresh',
                'Figures require the fresh production dataset')
        for row in data['points']:
            require(row['channel'] == channel, 'Wrong channel in fresh points')
            key = (channel, row['steps'], row['p'])
            require(key not in found, 'Duplicate production cell')
            require(row.get('complete') is True, 'Incomplete production cell')
            n, k = row['shots'], row['failures']
            require(isinstance(n, int) and isinstance(k, int) and 0 <= k <= n,
                    'Invalid binomial counts')
            require(k == row['full_prefix_failures'], 'Failure aliases disagree')
            require(n == expected.get(key), 'Production shot count differs from frozen plan')
            found[key] = n
        require('fits' in data and 'profiles' in data and 'failures' in data,
                'Analysis report is incomplete')
        settings = plan['analysis'][short]
        planned_fits = {f'Lmin{minimum}-{window}-{model}'
                        for minimum in settings['min_sizes'] for window in settings['windows']
                        for model in settings['models']}
        reported_fits = [r['fit_id'] for r in data['fits']]
        failed_fits = [r['fit_id'] for r in data['failures'] if r['stage'] == 'fit']
        require(len(reported_fits + failed_fits) == len(planned_fits)
                and set(reported_fits + failed_fits) == planned_fits,
                'Every planned fit must be saved or explicitly recorded as failed')
    require(found == expected, 'Production grid differs from frozen plan')
    require(sum(found.values()) == plan['total_new_fault_records'], 'Wrong total record count')
    x = report['fresh']['x']
    anchors = x['primary_crossing']['contrasts']
    primary = plan['analysis']['x']['primary_crossing']
    expected_anchors = {(*primary['size_pair'], p) for p in primary['rates']}
    require({(r['L_small'], r['L_large'], r['p']) for r in anchors} == expected_anchors
            and len(anchors) == len(expected_anchors), 'Incomplete primary contrasts')
    sizes = plan['analysis']['x']['sizes']
    expected_bands = {(a, b, p) for i, a in enumerate(sizes) for b in sizes[i + 1:]
                      for p in plan['analysis']['x']['rates']}
    bands = x['secondary_crossing']['bands']
    require({(r['L_small'], r['L_large'], r['p']) for r in bands} == expected_bands
            and len(bands) == len(expected_bands), 'Incomplete secondary crossing bands')
    for row in [*anchors, *bands]:
        require(all(finite(row[k]) for k in ('difference', 'ci_low', 'ci_high'))
                and row['ci_low'] <= row['difference'] <= row['ci_high'],
                'Invalid saved contrast interval')
    # The analysis writes a completion timestamp after every planned analysis.
    z = report['fresh']['z']
    settings = plan['analysis']['z']
    planned_profiles = {f'Lmin{minimum}-{window}-{model}'
                        for minimum in settings['min_sizes'] for window in settings['windows']
                        for model in settings['profile_models']}
    profiles = [r['fit_id'] for r in z['profiles']]
    failed_profiles = [r['fit_id'] for r in z['failures'] if r['stage'] == 'profile']
    require(len(profiles + failed_profiles) == len(planned_profiles)
            and set(profiles + failed_profiles) == planned_profiles,
            'Every planned profile must be saved or explicitly recorded as failed')
    require(finite(report.get('completed_at_unix'))
            and report.get('status') in ('complete', 'complete_with_numerical_failures'),
            'Analysis has not completed')


def wilson_interval(k, n):
    """Pointwise two-sided 95% Wilson binomial interval; no fitted quantities."""
    z = 1.959963984540054
    rate = k / n
    denominator = 1 + z * z / n
    center = (rate + z * z / (2 * n)) / denominator
    half = z * math.sqrt(rate * (1 - rate) / n + z * z / (4 * n * n)) / denominator
    # Clamp roundoff at k=0 and k=n so error-bar lengths cannot become negative.
    return min(rate, max(0., center - half)), max(rate, min(1., center + half))


def no_grid(fig):
    for ax in fig.axes:
        style_axes(ax)
        require(not any(g.get_visible() for g in ax.get_xgridlines() + ax.get_ygridlines()),
                'Gridlines must be disabled')
    return fig


def data_figure(data):
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    sizes = sorted({r['steps'] for r in data['points']})
    for size in sizes:
        rows = sorted((r for r in data['points'] if r['steps'] == size), key=lambda r: r['p'])
        p = np.array([r['p'] for r in rows])
        y = np.array([r['failures'] / r['shots'] for r in rows])
        bounds = np.array([wilson_interval(r['failures'], r['shots']) for r in rows])
        ax.errorbar(p, y, yerr=np.array([y - bounds[:, 0], bounds[:, 1] - y]),
                    color=COLORS[size], marker=MARKERS[size], linestyle='-',
                    linewidth=1.3, elinewidth=.9, capsize=2, markersize=4.4,
                    label=fr'$L={size}$')
    ax.set(xlabel=r'Gate-output fault probability $p$',
           ylabel='State failure probability\n(either of two logical parities)',
           title='X-only faults' if data['channel'] == 'x_only' else 'Z-only faults')
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=2))
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    if data['channel'] == 'x_only':
        ax.set_xticks(sorted({r['p'] for r in data['points']}))
    else:
        ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.margins(x=.035, y=.08)
    lo, hi = ax.get_ylim()
    ax.set_ylim(max(0., lo), min(1., hi))
    ax.legend(ncol=2, loc='best', fontsize=9, columnspacing=1.1, handlelength=2.1)
    fig.subplots_adjust(left=.13, right=.975, bottom=.16, top=.89)
    return no_grid(fig)


def x_crossing_figure(x, plan):
    fig = plt.figure(figsize=(11.0, 9.0))
    gs = fig.add_gridspec(3, 3, height_ratios=(1., 1., 1.), hspace=.65, wspace=.32)
    primary = x['primary_crossing']
    rows = sorted(primary['contrasts'], key=lambda r: r['p'])
    ax = fig.add_subplot(gs[0, :])
    y = np.array([r['difference'] for r in rows])
    error = np.array([[r['difference'] - r['ci_low'] for r in rows],
                      [r['ci_high'] - r['difference'] for r in rows]])
    ax.errorbar([r['p'] for r in rows], y, yerr=error, fmt='o', color='#0072B2',
                capsize=5, markersize=6, elinewidth=1.5)
    ax.axhline(0, color='.35', lw=.9, ls='--')
    pair = plan['analysis']['x']['primary_crossing']['size_pair']
    ax.set(title=fr'Primary: $L={pair[1]}$ minus $L={pair[0]}$, two planned probabilities',
           xlabel=r'Gate-output fault probability $p$', ylabel=r'$\Delta P_{\rm fail}$',
           xticks=[r['p'] for r in rows])
    ax.set_xlim(*plan['analysis']['x']['windows']['full'])
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=1))
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    low = min(0, min(r['ci_low'] for r in rows))
    high = max(0, max(r['ci_high'] for r in rows))
    pad = max(.002, (high - low) * .15)
    ax.set_ylim(low - pad, high + pad)

    secondary = x['secondary_crossing']
    pairs = sorted({(r['L_small'], r['L_large']) for r in secondary['bands']})
    limit = max(max(abs(r['ci_low']), abs(r['ci_high'])) for r in secondary['bands']) * 1.12
    limit = max(limit, .01)
    for index, pair in enumerate(pairs):
        ax = fig.add_subplot(gs[1 + index // 3, index % 3])
        rows = sorted((r for r in secondary['bands']
                       if (r['L_small'], r['L_large']) == pair), key=lambda r: r['p'])
        p = [r['p'] for r in rows]
        ax.fill_between(p, [r['ci_low'] for r in rows], [r['ci_high'] for r in rows],
                        color='#0072B2', alpha=.16)
        ax.plot(p, [r['difference'] for r in rows], 'o-', color='#0072B2', ms=3.2, lw=1.2)
        ax.axhline(0, color='.35', lw=.8, ls='--')
        ax.set(title=fr'Secondary: $L={pair[1]}$ minus $L={pair[0]}$',
               ylim=(-limit, limit), xlabel=r'Fault probability $p$')
        if index % 3 == 0:
            ax.set_ylabel(r'$\Delta P_{\rm fail}$')
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.tick_params(labelsize=9)
        ax.margins(x=.03)
    components = secondary['secondary_all_pairs']['components']
    intervals = ', '.join(f'[{a:.2%}, {b:.2%}]' for a, b in components) or 'empty'
    fig.suptitle('X: fresh-sample crossing contrasts', fontsize=14, y=.98)
    fig.subplots_adjust(left=.085, right=.975, top=.91, bottom=.18)
    fig.text(.085, .112, 'Familywise 95% intervals: four planned cells (primary); '
             f"all {secondary['number_cells']} cells (secondary).", fontsize=9)
    region_note = 'All-pair zero-compatible region: ' + intervals + '. Secondary bounds use straight-line interpolation.'
    fig.text(.085, .084, textwrap.fill(region_note, width=130), fontsize=8.5, va='top')
    fig.text(.085, .018, 'Crossing compatibility does not establish an asymptotic threshold.', fontsize=8.5)
    return no_grid(fig)


def fit_lookup(z):
    result = {}
    for fit in z['fits']:
        key = (fit['min_size'], fit['window'], fit['model'])
        require(key not in result, 'Duplicate saved fit')
        result[key] = fit
    return result


def fit_flags(fit):
    poor = finite(fit.get('deviance_p_value')) and fit['deviance_p_value'] < .01
    flagged = bool(fit.get('bound_hits')) or fit.get('optimizer_success') is not True \
        or fit.get('monotone_increasing') is not True
    return poor, flagged


def window_label(name, plan):
    lo, hi = plan['analysis']['z']['windows'][name]
    return f'{name.title()} window: {lo:.2%} to {hi:.2%}'


def z_sensitivity_figure(z, plan):
    fits = fit_lookup(z)
    minima = plan['analysis']['z']['min_sizes']
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.4), sharey=True)
    for ax, window in zip(axes, ('full', 'narrow')):
        missing = []
        for model, (label, color, marker, linestyle) in MODELS.items():
            rows = [fits.get((minimum, window, model)) for minimum in minima]
            values = [r['pc'] if r is not None and finite(r.get('pc')) else np.nan for r in rows]
            ax.plot(minima, values, color=color, linestyle=linestyle, lw=1.4)
            for minimum, value, fit in zip(minima, values, rows):
                if not np.isfinite(value):
                    missing.append(label + fr', $L_{{\min}}={minimum}$')
                    continue
                poor, flagged = fit_flags(fit)
                ax.plot(minimum, value, marker=marker, ms=6, mec=color,
                        mfc='white' if poor else color, mew=1.2, ls='none')
                if flagged:
                    ax.plot(minimum, value, marker='x', color='black', ms=9, mew=1.5, ls='none')
        if missing:
            ax.text(.03, .03, 'Unavailable:\n' + '\n'.join(missing), transform=ax.transAxes,
                    fontsize=8, va='bottom')
        ax.set(xlabel=r'Minimum included size $L_{\min}$', xticks=minima,
               xlim=(min(minima) - .4, max(minima) + .4), title=window_label(window, plan))
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=2))
    axes[0].set_ylabel(r'Fitted critical probability $p_c$')
    handles = [Line2D([], [], color=color, marker=marker, linestyle=linestyle, label=label)
               for label, color, marker, linestyle in MODELS.values()]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53, .12), ncol=2)
    fig.suptitle('Z: size-cutoff and fit-window sensitivity', fontsize=14, y=.98)
    fig.subplots_adjust(left=.10, right=.975, bottom=.28, top=.85, wspace=.15)
    fig.text(.10, .069, 'Point estimates only. Open markers: deviance p < 0.01. '
             'Black x: bound, monotonicity, or optimizer flag.', fontsize=8.5)
    fig.text(.10, .029, 'Two planned models shown; all six model families remain in the analysis report.', fontsize=8.5)
    return no_grid(fig)


def z_profile_figure(z, plan):
    fits = fit_lookup(z)
    profiles = {}
    for row in z['profiles']:
        key = (row['min_size'], row['window'], row['model'])
        require(key not in profiles, 'Duplicate saved profile')
        profiles[key] = row
    minima = plan['analysis']['z']['min_sizes']
    fig, axes = plt.subplots(len(minima), 2, figsize=(10.8, 9.4), sharey=True, squeeze=False)
    for row_index, minimum in enumerate(minima):
        for column, window in enumerate(('full', 'narrow')):
            ax = axes[row_index, column]
            notes = []
            cutoffs = set()
            for model, (label, color, marker, linestyle) in MODELS.items():
                key = (minimum, window, model)
                profile, fit = profiles.get(key), fits.get(key)
                if (profile is None or profile.get('critical_helper_exception')
                        or not profile.get('evaluations') or not finite(profile.get('lr_cutoff'))):
                    notes.append(label + ': profile unavailable')
                    continue
                cutoff = profile['lr_cutoff']
                cutoffs.add(cutoff)
                evaluations = sorted(profile['evaluations'], key=lambda r: r['pc'])
                values = [max(0., r['lr']) if r['optimizer_success'] and finite(r.get('lr'))
                          else np.nan for r in evaluations]
                ax.plot([r['pc'] for r in evaluations], values, color=color,
                        linestyle=linestyle, lw=1.5)
                agreement = profile.get('minimum_agreement') is True
                poor, flagged = fit_flags(fit) if fit is not None else (False, True)
                ax.plot(profile['free_pc'], 0, marker=marker, ms=5, mec=color,
                        mfc='white' if poor else color, clip_on=False, zorder=4)
                if flagged or not agreement:
                    ax.plot(profile['free_pc'], .35, marker='x', color='black', ms=7, ls='none')
                if poor:
                    notes.append(label + ': poor fit')
                if not agreement:
                    notes.append(label + ': minimum mismatch')
                if agreement:
                    for region in profile['regions']:
                        ax.axvspan(region['lower'], region['upper'], color=color, alpha=.07, zorder=-1)
                        for side, mark in [('lower', '<'), ('upper', '>')]:
                            if region.get(side + '_bound_truncated'):
                                ax.plot(region[side], cutoff, marker=mark, color=color, ms=6,
                                        clip_on=False, ls='none')
                            if region.get(side + '_numerical_gap'):
                                ax.plot(region[side], cutoff, marker='x', color=color, ms=6, ls='none')
                failed_p = [r['pc'] for r in evaluations
                            if not r['optimizer_success'] or not finite(r.get('lr'))]
                if failed_p:
                    ax.plot(failed_p, [11.4] * len(failed_p), marker='x', color=color, ms=4, ls='none')
            for cutoff in cutoffs:
                ax.axhline(cutoff, color='.35', ls=':', lw=1.)
            if notes:
                ax.text(.025, .965, '\n'.join(notes), transform=ax.transAxes, va='top', fontsize=8,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=2))
            ax.set(title=fr'$L_{{\min}}={minimum}$; ' + window_label(window, plan),
                   xlabel=r'Fixed critical probability $p_c$', ylim=(0, 12),
                   xlim=tuple(plan['analysis']['z']['windows'][window]))
            ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=2))
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.tick_params(labelsize=10)
    handles = [Line2D([], [], color=color, linestyle=linestyle, label=label)
               for label, color, marker, linestyle in MODELS.values()]
    handles.append(Line2D([], [], color='.35', linestyle=':', label='Approximate 95% cutoff'))
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .075), ncol=3, fontsize=9)
    fig.suptitle('Z: conditional likelihood profiles from fresh samples', fontsize=14, y=.99)
    fig.supylabel(r'$2\Delta$ negative log-likelihood', x=.025, y=.55, fontsize=13)
    fig.subplots_adjust(left=.10, right=.95, top=.93, bottom=.17, hspace=.56, wspace=.20)
    fig.text(.10, .052, 'Shading: detected conditional 95% regions; other parameters reoptimized. '
             'Model uncertainty is excluded.', fontsize=8.5)
    fig.text(.10, .027, 'Triangles: bound-truncated regions. Colored x: numerical gaps. '
             'Black x: fit flag or minimum mismatch.', fontsize=8.5)
    return no_grid(fig)


def build_figures(report, plan):
    validate_report(report, plan)
    x, z = report['fresh']['x'], report['fresh']['z']
    with paper_theme():
        return dict(zip(OUTPUTS, (data_figure(x), data_figure(z), x_crossing_figure(x, plan),
                                  z_sensitivity_figure(z, plan), z_profile_figure(z, plan))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=HERE / 'analysis' / 'report.json')
    parser.add_argument('--selection', type=Path, default=HERE / 'selection.json')
    parser.add_argument('--output', type=Path, default=HERE / 'figures')
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    plan = json.loads(args.selection.read_text())
    figures = build_figures(report, plan)
    args.output.mkdir(parents=True, exist_ok=True)
    with paper_theme():
        for name, fig in figures.items():
            fig.savefig(args.output / name, format='pdf', metadata={
                'Title': name.removesuffix('.pdf'),
                'Subject': 'H+6 full-prefix targeted follow-up; fresh production data only. '
                           'State failure means either of the two logical parities fails. '
                           'Data-plot error bars are pointwise 95% Wilson intervals.',
            })
            plt.close(fig)
    print(f'Saved {len(figures)} no-grid PDFs to {args.output.resolve()}')


if __name__ == '__main__':
    main()
