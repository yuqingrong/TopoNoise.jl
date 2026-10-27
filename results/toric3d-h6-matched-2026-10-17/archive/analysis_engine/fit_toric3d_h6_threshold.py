#!/usr/bin/env python3
"""H+6-CNOT final-only X/Z inference with mandatory publication guards.

Only the numerical binomial likelihood implementation is shared with the earlier
feedback analysis. Input validation, provenance, grid selection, acceptance, and
figures are specific to this H+6 experiment. Numerical candidates are never
threshold estimates unless every listed guard passes.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import struct
import sys
import tomllib

import numpy as np
from scipy.special import expit

ENGINE_PATH = Path(__file__).with_name('fit_toric3d_diagonal_feedback.py')
SPEC = importlib.util.spec_from_file_location('_h6_likelihood_engine', ENGINE_PATH)
ENGINE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ENGINE)
CHANNELS = ('x_only', 'z_only')
DECODERS = {'x_only': 'BP+OSD', 'z_only': 'unit_data_MWPM'}
PROTOCOL = dict(protocol='legacy_H6_final_only', circuit='plaquette-transfer H+6-CNOT',
                decoding='final-only', intermediate_feedback=False, feedback=False,
                final_decode_count=1, geometry='T=L')
NU_BOUNDS = (.3, 5.)
GUARD_POLICY = {
    'minimum_sizes': 4, 'minimum_rates': 5, 'goodness_alpha': .01,
    'minimum_bootstrap_replicates': 20, 'bootstrap_acceptance_fraction': .8,
    'crossing_cluster': 'exactly one connected common intersection of sign-changing brackets across every adjacent-size pair',
    'pc_sensitivity_tolerance': 'max(10% of primary rate window, bootstrap 95% interval width)',
    'nu_sensitivity_tolerance': 'max(25% of primary nu, bootstrap 95% interval width)',
    'sensitivity_variants': ['drop_smallest', 'narrow_window', 'cubic', 'inverse_L_correction'],
    'correction': 'add c*L_reference/L to the logit; fixed irrelevant exponent omega=1',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson_interval(failures, shots):
    low, high = ENGINE.wilson_interval(failures, shots)
    # Exact boundary values avoid tiny signed roundoff in plotted error bars.
    return (0.0 if failures == 0 else low, 1.0 if failures == shots else high)


def point_seed(seed, channel, size, p):
    bits = struct.pack('>d', p).hex()
    return int(hashlib.sha256(f'{seed}|{channel}|{size}|{bits}'.encode()).hexdigest()[:15], 16)


def integer(value, name):
    try:
        if isinstance(value, str):
            # Preserve exact 60-bit seeds rather than round-trip through float.
            result = int(value)
        else:
            result = int(value)
            if result != value:
                raise ValueError()
    except (ValueError, OverflowError, TypeError) as error:
        raise ValueError(f'{name} must be an integer') from error
    return result


def validate_rows(rows, *, require_both=True):
    """Require both channels, exact rectangular cubic grids, and valid unions."""
    output, seen = [], set()
    for original in rows:
        row = dict(original)
        channel = row.get('channel')
        if channel not in CHANNELS:
            raise ValueError('channel must be x_only or z_only')
        for field in ('steps', 'Ly', 'Lz', 'shots', 'failures',
                      'logical_1_failures', 'logical_2_failures', 'seed'):
            row[field] = integer(row.get(field), field)
        if not row['steps'] == row['Ly'] == row['Lz'] or row['steps'] < 2:
            raise ValueError('H+6 threshold geometry requires T=L and Ly=Lz=L >= 2')
        row['p'] = float(row['p'])
        if not math.isfinite(row['p']) or not 0 < row['p'] < .5:
            raise ValueError('p must be finite and in (0,.5)')
        n, k = row['shots'], row['failures']
        a, b = row['logical_1_failures'], row['logical_2_failures']
        if n <= 0 or not 0 <= k <= n:
            raise ValueError('counts require 0 <= failures <= shots and shots > 0')
        if not 0 <= a <= k or not 0 <= b <= k or k > a + b:
            raise ValueError('logical counts must describe a valid two-logical union')
        if str(row.get('complete')).lower() != 'true':
            raise ValueError('all fitted points must be complete')
        if row.get('decoder') != DECODERS[channel]:
            raise ValueError(f'unexpected decoder for {channel}')
        elapsed = float(row['elapsed_seconds'])
        if not math.isfinite(elapsed) or elapsed < 0:
            raise ValueError('elapsed_seconds must be finite and nonnegative')
        key = channel, row['steps'], row['p']
        if key in seen:
            raise ValueError(f'duplicate point {key}')
        seen.add(key)
        row['none_failures'] = k  # numerical-engine adapter; no feedback parser used
        output.append(row)
    if not output or (require_both and {r['channel'] for r in output} != set(CHANNELS)):
        raise ValueError('study must contain both X and Z channels')
    for channel in CHANNELS:
        selected = [r for r in output if r['channel'] == channel]
        sizes = {r['steps'] for r in selected}
        rates = {r['p'] for r in selected}
        if {(channel, L, p) for L in sizes for p in rates} != {k for k in seen if k[0] == channel}:
            raise ValueError(f'{channel} grid must be complete and rectangular')
    return sorted(output, key=lambda r: (r['channel'], r['steps'], r['p']))


def load_study(folder):
    folder = Path(folder)
    raw, metadata_path = folder / 'raw.csv', folder / 'metadata.toml'
    metadata = tomllib.loads(metadata_path.read_text())
    if metadata.get('raw_csv_sha256') != sha256(raw):
        raise ValueError('raw CSV SHA256 hash does not match metadata')
    if any(metadata.get(k) != v for k, v in PROTOCOL.items()):
        raise ValueError('protocol must be H+6 final-only, no intermediate feedback, T=L')
    if metadata.get('complete') is not True:
        raise ValueError('study metadata must declare complete=true')
    with raw.open(newline='') as stream:
        rows = validate_rows(list(csv.DictReader(stream)), require_both=False)
    config = metadata.get('config', {})
    if not config.get('channels') or not set(config['channels']).issubset(CHANNELS):
        raise ValueError('configuration must specify valid channels')
    expected = set()
    for channel in config['channels']:
        settings = config['channels'][channel]
        sizes, rates = settings['sizes'], settings['rates']
        if len(sizes) != len(set(sizes)) or len(rates) != len(set(rates)):
            raise ValueError('configuration grid contains duplicates')
        expected.update((channel, L, p) for L in sizes for p in rates)
    if expected != {(r['channel'], r['steps'], r['p']) for r in rows}:
        raise ValueError('raw grid does not exactly match configured grid')
    if any(r['shots'] != config.get('shots') for r in rows):
        raise ValueError('raw shot counts do not match configuration')
    for row in rows:
        if row['seed'] != point_seed(config['seed'], row['channel'], row['steps'], row['p']):
            raise ValueError('raw point seed does not match deterministic assignment')
    if metadata.get('expected_points') != len(rows) or metadata.get('completed_points') != len(rows):
        raise ValueError('metadata point counts disagree with raw grid')
    source_hashes = metadata.get('source_hashes', {})
    if not source_hashes:
        raise ValueError('source_hashes provenance is required')
    root = Path(__file__).resolve().parents[1]
    for path, digest in source_hashes.items():
        resolved = (root / path).resolve()
        if not resolved.is_relative_to(root) or not resolved.is_file() or sha256(resolved) != digest:
            raise ValueError(f'source hash mismatch: {path}')
    config_path = folder / 'config.toml'
    if not config_path.is_file() or metadata.get('config_sha256') != sha256(config_path):
        raise ValueError('config SHA256 hash does not match saved config.toml')
    saved_config = tomllib.loads(config_path.read_text())
    if saved_config != config:
        raise ValueError('saved config.toml differs from metadata config')
    return rows, metadata


def load_studies(folders):
    """Combine disjoint channel studies only; each is validated before pooling."""
    rows, metadata, channels = [], [], set()
    for folder in folders:
        new_rows, information = load_study(folder)
        new_channels = {row['channel'] for row in new_rows}
        if channels & new_channels:
            raise ValueError('input studies have channel overlap; duplicate/pooling not permitted')
        channels.update(new_channels)
        rows.extend(new_rows)
        metadata.append(information)
    for key in ('source_hashes', 'backend_versions', 'decoder_settings'):
        present = [information.get(key) for information in metadata]
        if any(value != present[0] for value in present[1:]):
            raise ValueError(f'input studies have incompatible {key}')
    return validate_rows(rows), metadata


def crossing_clusters(points):
    """Keep all empirical crossings and require common bracketing evidence."""
    sizes = sorted({r['steps'] for r in points})
    rates = sorted({r['p'] for r in points})
    by_key = {(r['steps'], r['p']): r['failures'] / r['shots'] for r in points}
    brackets, all_crossings = [], []
    for a, b in zip(sizes[:-1], sizes[1:]):
        differences = [by_key[b, p] - by_key[a, p] for p in rates]
        pair = []
        for i in range(len(rates) - 1):
            left, right = differences[i:i+2]
            if left * right < 0:
                crossing = rates[i] - left * (rates[i+1] - rates[i]) / (right - left)
                pair.append(dict(sizes=[a, b], p=crossing, bracket=[rates[i], rates[i+1]]))
        for i in range(1, len(rates) - 1):
            if differences[i] == 0 and differences[i-1] * differences[i+1] < 0:
                pair.append(dict(sizes=[a, b], p=rates[i], bracket=[rates[i-1], rates[i+1]]))
        brackets.append(pair)
        all_crossings.extend(pair)
    common = [[rates[0], rates[-1]]] if rates else []
    for pair in brackets:
        common = [[max(a, item['bracket'][0]), min(b, item['bracket'][1])]
                  for a, b in common for item in pair
                  if max(a, item['bracket'][0]) <= min(b, item['bracket'][1])]
    merged = []
    for low, high in sorted(set(map(tuple, common))):
        if merged and low <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    return all_crossings, merged


def candidate_fit(points, degree=2, correction=False, warm_start=None):
    rates = sorted({r['p'] for r in points})
    if len(rates) < 2 or len(points) <= degree + 3 + int(correction):
        return None
    try:
        return ENGINE._candidate_fit(
            points, points[0]['channel'], 'none', p_min=rates[0], p_max=rates[-1],
            degree=degree, nu_bounds=NU_BOUNDS, goodness_alpha=.01,
            correction=correction, warm_start=warm_start)
    except (ArithmeticError, RuntimeError, ValueError):
        return None


def numerical_guards(candidate):
    if candidate is None:
        return dict(interior_optimizer=False, monotonicity=False, goodness_of_fit=False)
    return dict(
        interior_optimizer=bool(candidate['optimizer_success'] and not candidate['bound_hits']
                                and candidate['p_min'] < candidate['pc'] < candidate['p_max']),
        monotonicity=bool(candidate['monotone_increasing']),
        goodness_of_fit=bool(candidate['degrees_of_freedom'] > 0
                             and candidate['deviance_p_value'] >= .01),
    )


def bootstrap_candidate(points, candidate, replicates, seed):
    samples = []
    rng = np.random.default_rng(seed)
    if candidate is not None:
        for index in range(replicates):
            sampled = []
            for row, probability in zip(points, candidate['fitted_probabilities']):
                k = int(rng.binomial(row['shots'], probability))
                sampled.append(dict(row, failures=k, none_failures=k))
            fitted = candidate_fit(sampled, warm_start=candidate['parameters'])
            guards = numerical_guards(fitted)
            usable = guards['interior_optimizer'] and guards['monotonicity']
            samples.append(dict(replicate=index, accepted=usable,
                                pc=fitted['pc'] if fitted is not None else None,
                                nu=fitted['nu'] if fitted is not None else None,
                                failed_guards=[k for k, v in guards.items()
                                               if not v and k != 'goodness_of_fit']))
    accepted = [r for r in samples if r['accepted']]
    intervals = {key: np.quantile([r[key] for r in accepted], [.025, .975]).tolist()
                 if len(accepted) >= 2 else None for key in ('pc', 'nu')}
    return dict(replicates=replicates, accepted=len(accepted), seed=seed,
                pc_interval=intervals['pc'], nu_interval=intervals['nu'], samples=samples,
                method='seeded parametric binomial bootstrap conditional on candidate model; numerical acceptance only')


def analyze_channel(rows, channel, *, bootstrap=200, seed=20261015):
    if bootstrap < 0:
        raise ValueError('bootstrap must be nonnegative')
    rows = validate_rows(rows)
    points = [r for r in rows if r['channel'] == channel]
    if not points:
        raise ValueError('unknown channel')
    sizes, rates = sorted({r['steps'] for r in points}), sorted({r['p'] for r in points})
    candidate = candidate_fit(points)
    crossings, clusters = crossing_clusters(points)
    guards = dict(grid_resolution=len(sizes) >= 4 and len(rates) >= 5)
    guards.update(numerical_guards(candidate))
    guards['crossings'] = bool(len(sizes) >= 2 and len(clusters) == 1 and candidate is not None
                               and any(a <= candidate['pc'] <= b for a, b in clusters))
    # A failed candidate is still bootstrapped diagnostically; samples cannot override guards.
    boot = bootstrap_candidate(points, candidate, bootstrap, seed)
    guards['bootstrap'] = bool(
        bootstrap >= GUARD_POLICY['minimum_bootstrap_replicates']
        and boot['accepted'] >= math.ceil(.8 * bootstrap)
        and boot['pc_interval'] is not None
        and rates[0] < boot['pc_interval'][0] < boot['pc_interval'][1] < rates[-1]
        and NU_BOUNDS[0] < boot['nu_interval'][0] < boot['nu_interval'][1] < NU_BOUNDS[1])
    variants = [
        ('drop_smallest', [r for r in points if r['steps'] != sizes[0]], 2, False),
        ('narrow_window', [r for r in points if rates[0] < r['p'] < rates[-1]], 2, False),
        ('cubic', points, 3, False),
        ('inverse_L_correction', points, 2, True),
    ]
    sensitivities = []
    for name, selected, degree, correction in variants:
        enough = len({r['steps'] for r in selected}) >= 4 and len({r['p'] for r in selected}) >= 5
        fitted = candidate_fit(selected, degree, correction) if enough else None
        checks = numerical_guards(fitted)
        stable = False
        pc_tolerance = nu_tolerance = None
        if fitted is not None and candidate is not None:
            pc_tolerance = max(.1 * (rates[-1] - rates[0]),
                               np.ptp(boot['pc_interval']) if boot['pc_interval'] else 0.)
            nu_tolerance = max(.25 * candidate['nu'],
                               np.ptp(boot['nu_interval']) if boot['nu_interval'] else 0.)
            stable = (abs(fitted['pc'] - candidate['pc']) <= pc_tolerance
                      and abs(fitted['nu'] - candidate['nu']) <= nu_tolerance)
        passed = bool(enough and all(checks.values()) and stable)
        guards[f'sensitivity_{name}'] = passed
        sensitivities.append(dict(name=name, passed=passed, candidate=fitted,
                                  sufficient_grid=enough, guards=checks, stable=bool(stable),
                                  pc_tolerance=pc_tolerance, nu_tolerance=nu_tolerance))
    accepted = all(guards.values())
    return dict(channel=channel, status='accepted' if accepted else 'diagnostic_only',
                pc=candidate['pc'] if accepted else None,
                nu=candidate['nu'] if accepted else None,
                candidate=candidate, guards=guards,
                failed_guards=[key for key, passed in guards.items() if not passed],
                sizes=sizes, rates=rates, crossings=crossings,
                crossing_common_brackets=clusters, bootstrap=boot, sensitivity=sensitivities)


def make_figure(rows, fits, output):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/toponoise-h6-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    figure, axes = plt.subplots(1, 2, figsize=(13.6, 7.1))
    for axis, channel in zip(axes, CHANNELS):
        fit = fits[channel]
        selected = [r for r in rows if r['channel'] == channel]
        colors = plt.cm.viridis(np.linspace(.1, .85, len(fit['sizes'])))
        for color, size in zip(colors, fit['sizes']):
            points = sorted([r for r in selected if r['steps'] == size], key=lambda r: r['p'])
            observed = np.array([r['failures'] / r['shots'] for r in points])
            intervals = np.array([wilson_interval(r['failures'], r['shots']) for r in points])
            axis.errorbar([100 * r['p'] for r in points], observed,
                          yerr=np.maximum(0, [observed - intervals[:, 0], intervals[:, 1] - observed]),
                          marker='o', ms=4, lw=1.1, capsize=2, color=color, label=f'L={size}')
        axis.set(xlabel='Gate-output error probability p (%)', ylabel='Either-logical failure probability',
                 title=f'{channel[0].upper()} channel · ' + ('BP+OSD' if channel == 'x_only' else 'MWPM'))
        maximum = max(wilson_interval(r['failures'], r['shots'])[1] for r in selected)
        axis.set_ylim(0, min(1, 1.28 * maximum))
        axis.set_title(axis.get_title(), pad=40)
        axis.grid(alpha=.15)
        axis.legend(loc='lower center', bbox_to_anchor=(.5, 1.005), ncol=3,
                    fontsize=8, frameon=False, borderaxespad=0)
        inset = axis.inset_axes([.065, .49, .46, .40])
        candidate = fit['candidate']
        accepted = fit['status'] == 'accepted'
        inset.set_title('Accepted scaling fit' if accepted else 'Diagnostic scaling fit', fontsize=9)
        if candidate is not None:
            all_x = []
            for color, size in zip(colors, fit['sizes']):
                points = [r for r in selected if r['steps'] == size]
                x = [(r['p'] - candidate['pc']) / candidate['p_scale']
                     * (size / candidate['L_reference']) ** (1 / candidate['nu']) for r in points]
                all_x.extend(x)
                inset.plot(x, [r['failures']/r['shots'] for r in points], 'o', ms=2.5, color=color)
            grid = np.linspace(min(all_x), max(all_x), 300)
            inset.plot(grid, expit(np.polynomial.polynomial.polyval(grid, candidate['coefficients'])),
                       color='#222222', lw=1.5)
            inset.set_xlabel('Scaled p − candidate p₀', fontsize=8, labelpad=2)
            inset.set_ylabel('Failure', fontsize=8, labelpad=1)
            inset.tick_params(labelsize=7)
            inset.grid(alpha=.15)
        else:
            inset.text(.5, .5, 'Fit unavailable', ha='center', va='center', transform=inset.transAxes)
        if accepted:
            interval = fit['bootstrap']['pc_interval']
            label = f'Threshold p_c = {100*fit["pc"]:.4f}% (95% CI {100*interval[0]:.4f}–{100*interval[1]:.4f}%)'
        else:
            label = 'Threshold withheld · ' + ', '.join(fit['failed_guards'][:2]).replace('_', ' ')
            if len(fit['failed_guards']) > 2:
                label += f' (+{len(fit["failed_guards"])-2} guards)'
        axis.text(.02, -.18, label, transform=axis.transAxes, fontsize=9,
                  color='#137555' if accepted else '#a23b32')
    figure.suptitle('Legacy H+6-CNOT toric code · T=L · one final decode · no feedback', fontsize=14, y=.98)
    figure.text(.5, .022, '95% Wilson error bars. Insets: binomial-logit finite-size scaling candidates. '
                'X uses an ideal logical |00> seed; Z uses a noisy-H logical |++> seed.',
                ha='center', fontsize=8.5)
    figure.subplots_adjust(left=.065, right=.985, top=.81, bottom=.25, wspace=.24)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf', 'svg'):
        figure.savefig(output / f'h6-threshold-failure-rates.{extension}', dpi=200)
    return figure


def write_csv(path, rows, fieldnames):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row[key]) if isinstance(row.get(key), (dict, list))
                             else row.get(key) for key in fieldnames})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, type=Path, action='append')
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--bootstrap', type=int, default=200)
    parser.add_argument('--seed', type=int, default=20261015)
    args = parser.parse_args(argv)
    rows, metadata = load_studies(args.input_dir)
    fits = {channel: analyze_channel(rows, channel, bootstrap=args.bootstrap, seed=args.seed + index)
            for index, channel in enumerate(CHANNELS)}
    versions = {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'matplotlib')}
    report = dict(format_version=1, experiment=PROTOCOL, input_metadata=metadata,
                  input_hashes=[dict(directory=str(folder.resolve()),
                                     **{name: sha256(folder / name)
                                        for name in ('raw.csv', 'metadata.toml', 'config.toml')})
                                for folder in args.input_dir],
                  analysis_hashes={str(Path(__file__).name): sha256(__file__),
                                   ENGINE_PATH.name: sha256(ENGINE_PATH)},
                  versions=dict(python=platform.python_version(), **versions),
                  python_executable=sys.executable,
                  thread_environment={name: os.environ.get(name) for name in
                                      ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
                  guard_policy=GUARD_POLICY,
                  bootstrap_seed=args.seed, fits=fits,
                  interpretation='Thresholds are decoder- and preparation-specific. Diagnostic candidates are not threshold estimates.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'threshold-report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    write_csv(args.output_dir / 'bootstrap-samples.csv',
              [dict(channel=channel, **sample) for channel, fit in fits.items() for sample in fit['bootstrap']['samples']],
              ['channel', 'replicate', 'accepted', 'pc', 'nu', 'failed_guards'])
    write_csv(args.output_dir / 'sensitivity-fits.csv',
              [dict(channel=channel, **sensitivity) for channel, fit in fits.items() for sensitivity in fit['sensitivity']],
              ['channel', 'name', 'passed', 'sufficient_grid', 'stable', 'pc_tolerance', 'nu_tolerance', 'guards', 'candidate'])
    write_csv(args.output_dir / 'threshold-points.csv',
              [dict(row, rate=row['failures']/row['shots'],
                    ci_low=wilson_interval(row['failures'], row['shots'])[0],
                    ci_high=wilson_interval(row['failures'], row['shots'])[1]) for row in rows],
              [key for key in rows[0] if key != 'none_failures'] + ['rate', 'ci_low', 'ci_high'])
    make_figure(rows, fits, args.output_dir)
    print(json.dumps({channel: {k: fit[k] for k in ('status', 'pc', 'nu', 'failed_guards')}
                      for channel, fit in fits.items()}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
