#!/usr/bin/env python3
"""Paired four-arm H+6 Z analysis; all inputs are measured production counts."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import tomllib

import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'examples'))
SPEC = importlib.util.spec_from_file_location('h6_fit_policy', ROOT / 'examples/fit_toric3d_h6_threshold.py')
H6 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(H6)
ARMS = ('final_only', 'seed_only', 'one_step', 'window_complete')
LABELS = ('Final-only recovery', 'Seed correction only', 'Causal one-step Z feedback', 'Complete-window Z feedback')
COLORS = {3: '#7F7F7F', 5: '#0072B2', 7: '#E69F00', 9: '#009E73', 11: '#CC79A7'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value):
    result = int(value)
    require(str(result) == str(value) or result == value, 'noninteger count')
    return result


def paired(joint, first, second):
    require(len(joint) == 16 and all(n >= 0 and int(n) == n for n in joint) and sum(joint) > 0,
            'joint counts must be 16 nonnegative integers with positive total')
    a, b = 1 << ARMS.index(first), 1 << ARMS.index(second)
    require(a != b, 'paired arms must differ')
    groups = dict(neither=0, first_only=0, second_only=0, both=0)
    for mask, count in enumerate(joint):
        key = 'both' if mask & a and mask & b else 'first_only' if mask & a else 'second_only' if mask & b else 'neither'
        groups[key] += count
    n = sum(joint)
    delta = (groups['second_only'] - groups['first_only']) / n
    def interval(k):
        lo, hi = H6.ENGINE.wilson_interval(k, n, z=2.241402727604947)
        return (0. if k == 0 else lo, 1. if k == n else hi)
    a_ci = interval(groups['first_only'])
    b_ci = interval(groups['second_only'])
    return dict(first=first, second=second, **groups, difference=delta,
                standard_error=math.sqrt(max(0, ((groups['first_only'] + groups['second_only']) / n - delta**2) / (n - 1))) if n > 1 else None,
                ci_low=max(-1., b_ci[0] - a_ci[1]), ci_high=min(1., b_ci[1] - a_ci[0]))


def validate_point(original):
    require(str(original.get('complete')).lower() == 'true', 'incomplete production point')
    row = dict(original)
    for key in ('steps', 'Ly', 'Lz', 'shots', 'target_shots', 'seed'):
        row[key] = integer(row[key])
    require(row['steps'] == row['Ly'] == row['Lz'] >= 2, 'geometry must be T=L')
    require(row['shots'] == row['target_shots'] and row['shots'] > 0, 'incorrect shot count')
    row['p'], row['elapsed_seconds'] = float(row['p']), float(row['elapsed_seconds'])
    require(math.isfinite(row['p']) and 0 <= row['p'] < .5, 'invalid probability')
    require(math.isfinite(row['elapsed_seconds']) and row['elapsed_seconds'] >= 0, 'invalid timing')
    require(re.fullmatch('[0-9a-f]{16}', str(row['point_seed'])) is not None, 'seed must be exactly 16 hex digits')
    joint = [integer(row[f'joint_{mask:02d}']) for mask in range(16)]
    require(min(joint) >= 0 and sum(joint) == row['shots'], 'joint histogram does not sum to shots')
    for i, arm in enumerate(ARMS):
        keys = [arm + suffix for suffix in ('_logical_1_failures', '_logical_2_failures', '_failures')]
        a, b, union = [integer(row[k]) for k in keys]
        row.update(zip(keys, (a, b, union)))
        require(0 <= a <= union and 0 <= b <= union <= row['shots'] and union <= a + b, 'invalid logical union')
        require(union == sum(n for mask, n in enumerate(joint) if mask & (1 << i)), 'marginal differs from joint histogram')
    row['joint'] = joint
    return row


def load_points(directories):
    import importlib.metadata

    expected_protocol = dict(
        circuit='legacy H+6-CNOT plaquette transfer', channel='z_only',
        geometry='T=L; Ly=Lz=L', initial_state='logical |++>', arm_order=list(ARMS),
        final_only='no seed or intermediate correction; one final MWPM',
        seed_only='complete seed correction; no intermediate correction; one final MWPM',
        one_step='complete seed correction; causal local one-step intermediate correction; one final MWPM',
        window_complete='complete seed correction and full local-window recovery; one final MWPM',
        noise='independent Bernoulli-p Z faults on 3 L^2 seed-H outputs and 13 L^3 transfer outputs',
        sampler='sample_toric3d_diagonal_z_faults: seed bits first, then participating transfer outputs in schedule order',
        shared_final_decoder='family.prefixes[end], unit-weight data-edge PyMatching MWPM, syndrome-validated recovery',
        failure='either residual logical-Z parity after final recovery',
        joint_histogram_encoding='mask = final_only + 2*seed_only + 4*one_step + 8*window_complete; each term is its state-failure Boolean',
        point_order='sizes outer loop, rates inner loop; both strictly increasing',
        point_seed_generation='MersenneTwister(integer base seed); repeated scalar rand(rng,UInt64); discard duplicate values; 16 lowercase hexadecimal digits',
        point_rng='Random.Xoshiro(point_seed)',
        resume='recreate point Xoshiro and sample/discard completed records without decoding; rebuild raw CSV from authoritative checkpoint')
    sources = [ROOT / 'Project.toml', HERE / 'run_study.jl', *sorted((ROOT / 'src').rglob('*.jl'))]
    if (ROOT / 'Manifest.toml').exists():
        sources.append(ROOT / 'Manifest.toml')
    expected_sources = {str(path.relative_to(ROOT)): H6.sha256(path) for path in sources}
    expected_packages = {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'pymatching')}
    points, inputs = [], []
    shared_runtime = None
    for directory in map(Path, directories):
        require(not (directory / 'failed-record.toml').exists(), 'input has a failed simulation record')
        with (directory / 'raw.csv').open() as stream:
            rows = [validate_point(r) for r in csv.DictReader(stream)]
        checkpoint = tomllib.loads((directory / 'checkpoint.toml').read_text())
        metadata = tomllib.loads((directory / 'metadata.toml').read_text())
        require(rows and checkpoint.get('format_version') == metadata.get('format_version') == 1,
                'empty input or unsupported metadata/checkpoint format')
        for key in ('config', 'runtime', 'source_hashes', 'protocol', 'point_seeds_hex',
                    'generation', 'started_at_unix', 'updated_at_unix'):
            require(key in checkpoint and metadata.get(key) == checkpoint[key],
                    f'metadata/checkpoint mismatch: {key}')
        require(checkpoint['protocol'] == expected_protocol, 'unexpected H+6 Z four-arm protocol')
        require(checkpoint['source_hashes'] == expected_sources, 'simulation source hashes differ from current files')
        runtime = checkpoint['runtime']
        require(runtime.get('julia_threads') == runtime.get('blas_threads') == 1,
                'simulation must use one Julia and BLAS thread')
        environment = runtime.get('environment', {})
        require(all(environment.get(key) == '1' for key in
                    ('JULIA_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'))
                and environment.get('JULIA_CONDAPKG_BACKEND') == 'Null', 'unexpected simulation thread/runtime settings')
        require(Path(runtime.get('python_executable', '')).resolve() == Path(sys.executable).resolve()
                and Path(environment.get('JULIA_PYTHONCALL_EXE', '')).resolve() == Path(sys.executable).resolve()
                # Embedded libpython and the executable can have distinct build
                # timestamps in sys.version within the same Conda installation.
                and isinstance(runtime.get('python_version'), str)
                and runtime['python_version'].split()[:1] == sys.version.split()[:1]
                and runtime.get('packages') == expected_packages,
                'simulation Python executable/version/packages differ from analysis runtime')
        require(isinstance(runtime.get('julia_version'), str) and runtime['julia_version'], 'missing Julia version')
        require(shared_runtime is None or runtime == shared_runtime, 'simulation runtimes differ across shards')
        shared_runtime = runtime
        config = checkpoint['config']
        sizes, rates = config.get('sizes', []), config.get('rates', [])
        require(sizes and all(type(L) is int and L >= 2 for L in sizes)
                and sizes == sorted(set(sizes)), 'invalid configured sizes')
        require(rates and all(type(p) in (float, int) and math.isfinite(p) and 0 <= p < .5 for p in rates)
                and rates == sorted(set(rates)), 'invalid configured rates')
        require(config.get('dimensions') == [[L, L, L] for L in sizes], 'configured geometry is not T=L')
        require(all(type(config.get(key)) is int and config[key] > 0 for key in ('shots', 'batch_size'))
                and type(config.get('seed')) is int and config['seed'] >= 0, 'invalid configured shots/batch size/base seed')
        grid = [(L, p) for L in sizes for p in rates]
        states = checkpoint['states']
        require(len(states) == len(rows) == len(checkpoint['point_seeds_hex']) == len(grid),
                'checkpoint/raw/grid/seed point counts differ')
        require([(row['steps'], row['p']) for row in rows] == grid, 'raw row identity/order differs from configured grid')
        require(all(row['seed'] == config['seed'] and row['shots'] == row['target_shots'] == config['shots'] for row in rows),
                'raw base seed or shot target differs from config')
        require(metadata.get('complete') is True and metadata.get('expected_points') == metadata.get('completed_points') == len(grid)
                and metadata.get('completed_shots') == len(grid) * config['shots'], 'metadata does not describe a complete grid')
        require(metadata.get('raw_sha256') == H6.sha256(directory / 'raw.csv'), 'raw CSV hash mismatch')
        require(metadata.get('checkpoint_sha256') == H6.sha256(directory / 'checkpoint.toml'), 'checkpoint hash mismatch')
        generation, payload = checkpoint['generation'], checkpoint.get('payload_sha256', '')
        require(type(generation) is int and generation > 0 and re.fullmatch('[0-9a-f]{64}', payload) is not None,
                'invalid checkpoint generation/payload digest')
        require(all(type(checkpoint[key]) in (int, float) and math.isfinite(checkpoint[key]) for key in
                    ('started_at_unix', 'updated_at_unix')) and 0 < checkpoint['started_at_unix'] <= checkpoint['updated_at_unix'],
                'invalid checkpoint timestamps')
        snapshot_name = f'checkpoints/checkpoint-{generation:06d}-{payload[:12]}.toml'
        snapshot = directory / snapshot_name
        require(metadata.get('current_snapshot') == snapshot_name and snapshot.is_file()
                and snapshot.read_bytes() == (directory / 'checkpoint.toml').read_bytes(),
                'immutable checkpoint snapshot mismatch')
        # Julia owns its canonical TOML payload checksum and exact MersenneTwister
        # seed regeneration. This loader binds their recorded values to every row.
        for row, state, seed in zip(rows, states, checkpoint['point_seeds_hex']):
            require(type(state.get('completed_shots')) is int, 'checkpoint shot count must be an integer')
            for key, length in (('joint_counts', 16), ('failures', 4), ('logical_1_counts', 4), ('logical_2_counts', 4)):
                require(isinstance(state.get(key), list) and len(state[key]) == length
                        and all(type(n) is int and 0 <= n <= row['shots'] for n in state[key]),
                        f'invalid checkpoint count array: {key}')
            require(row['point_seed'] == seed and row['shots'] == state['completed_shots'], 'checkpoint seed/shots differ')
            require(row['elapsed_seconds'] == state.get('elapsed_seconds'), 'checkpoint elapsed time differs')
            require(row['joint'] == state['joint_counts'], 'checkpoint joint histogram differs')
            for i, arm in enumerate(ARMS):
                require(row[arm + '_failures'] == state['failures'][i], 'checkpoint marginal differs')
                for logical in (1, 2):
                    require(row[f'{arm}_logical_{logical}_failures'] == state[f'logical_{logical}_counts'][i], 'checkpoint logical count differs')
        points.extend(rows)
        inputs.append(dict(directory=str(directory.resolve()), metadata=metadata,
                           hashes={name: H6.sha256(directory / name)
                                   for name in ('raw.csv', 'checkpoint.toml', 'metadata.toml', snapshot_name)}))
    require(points, 'no input points')
    keys = [(r['steps'], r['p']) for r in points]
    require(len(set(keys)) == len(keys), 'duplicate production point')
    require(len({r['point_seed'] for r in points}) == len(points), 'duplicate point seed')
    sizes, rates = sorted({r['steps'] for r in points}), sorted({r['p'] for r in points})
    require(set(keys) == {(L, p) for L in sizes for p in rates}, 'production grid must be rectangular')
    return sorted(points, key=lambda r: (r['steps'], r['p'])), inputs


def adapted(points, arm):
    return [dict(channel='z_only', steps=r['steps'], p=r['p'], shots=r['shots'],
                 failures=r[arm + '_failures'], none_failures=r[arm + '_failures']) for r in points]


def validate_selection(points, inputs, selection, selection_path, *, bootstrap, seed):
    """Bind production inputs and inference settings to the recorded freeze."""
    study = Path(selection_path).resolve().parent
    require(selection.get('stage') == 'production' and selection.get('frozen_before_sampling') is True
            and selection.get('pilot_excluded') is True, 'production selection must be frozen and exclude the pilot')
    frozen = selection.get('frozen_at_unix')
    require(type(frozen) in (int, float) and math.isfinite(frozen) and frozen > 0,
            'invalid selection freeze timestamp')
    require(bootstrap == selection.get('bootstrap_replicates') and seed == selection.get('bootstrap_seed'),
            'bootstrap replicates/seed differ from frozen selection')
    require(type(bootstrap) is int and bootstrap >= H6.GUARD_POLICY['minimum_bootstrap_replicates']
            and type(seed) is int and seed >= 0, 'invalid frozen bootstrap settings')
    require(sorted({r['steps'] for r in points}) == selection.get('sizes'), 'sizes differ from frozen selection')
    require(sorted({r['p'] for r in points}) == selection.get('rates'), 'rates differ from frozen selection')
    require(all(r['shots'] == selection.get('shots') for r in points), 'shots differ from frozen selection')
    fit_rates = selection.get('fit_rates', [])
    require(fit_rates and all(type(p) in (int, float) and math.isfinite(p) and p > 0 for p in fit_rates)
            and fit_rates == sorted(set(fit_rates)) and set(fit_rates) <= set(selection['rates']),
            'primary fit rates must be an ordered positive subset of the frozen sampled rates')
    records = sum(row['shots'] for row in points)
    require(selection.get('total_fault_records') == records
            and selection.get('total_protocol_outcomes') == len(ARMS) * records,
            'sample totals differ from frozen selection')
    shards = selection.get('shards', {})
    require(isinstance(shards, dict) and shards and all(isinstance(name, str) and Path(name).name == name
            and name not in ('.', '..', 'pilot') for name in shards), 'invalid frozen production shards')
    actual = {Path(item['directory']).resolve(): item for item in inputs}
    require(len(actual) == len(inputs) and set(actual) == {(study / name).resolve() for name in shards},
            'input directories differ from frozen production shards')
    for name, declared in shards.items():
        metadata = actual[(study / name).resolve()]['metadata']
        config = metadata['config']
        require(config['sizes'] == declared.get('sizes') and config['seed'] == declared.get('seed')
                and config['rates'] == selection['rates'] and config['shots'] == selection['shots'],
                f'production shard configuration/master seed differs: {name}')
        require(frozen < metadata['started_at_unix'], 'selection was not frozen before production started')
        require(metadata['source_hashes'] == selection.get('pilot_source_sha256'),
                'production sources differ from the frozen pilot sources')
    numerical_files = ('examples/fit_toric3d_h6_threshold.py', 'examples/fit_toric3d_diagonal_feedback.py')
    numerical_hashes = {name: H6.sha256(ROOT / name) for name in numerical_files}
    require(selection.get('numerical_fit_source_sha256') == numerical_hashes,
            'numerical fit sources differ from frozen selection')
    pilot = study / 'pilot'
    pilot_names = ('raw.csv', 'metadata.toml', 'checkpoint.toml')
    require(all((pilot / name).is_file() for name in pilot_names), 'frozen pilot evidence is missing')
    pilot_hashes = {name: H6.sha256(pilot / name) for name in pilot_names}
    require(selection.get('pilot_input_sha256') == pilot_hashes, 'pilot input hashes differ from frozen selection')
    pilot_points, pilot_inputs = load_points([pilot])
    pilot_metadata = pilot_inputs[0]['metadata']
    require(pilot_metadata['source_hashes'] == selection.get('pilot_source_sha256'),
            'pilot source hashes differ from frozen selection')
    require(pilot_metadata['updated_at_unix'] <= frozen, 'production selection predates completed pilot evidence')
    require(not ({row['point_seed'] for row in pilot_points} & {row['point_seed'] for row in points}),
            'pilot fault-record streams overlap production')
    return dict(validated=True, selection_sha256=H6.sha256(selection_path),
                frozen_at_unix=frozen, bootstrap_replicates=bootstrap, bootstrap_seed=seed,
                pilot_input_sha256=pilot_hashes, numerical_fit_source_sha256=numerical_hashes)


def scaling_diagnostics(points, rates, bootstrap, seed):
    selected = [r for r in adapted(points, 'one_step') if r['p'] in rates and r['p'] > 0]
    sizes = sorted({r['steps'] for r in selected})
    rates = sorted({r['p'] for r in selected})
    require(len(sizes) >= 2 and len(rates) >= 3, 'too few points for scaling diagnostics')
    candidate = H6.candidate_fit(selected)
    crossings, clusters = H6.crossing_clusters(selected)
    guards = dict(grid_resolution=len(sizes) >= 4 and len(rates) >= 5, **H6.numerical_guards(candidate))
    guards['crossings'] = bool(len(clusters) == 1 and candidate is not None and
                              any(a <= candidate['pc'] <= b for a, b in clusters))
    boot = H6.bootstrap_candidate(selected, candidate, bootstrap, seed)
    guards['bootstrap'] = bool(bootstrap >= H6.GUARD_POLICY['minimum_bootstrap_replicates'] and
        boot['accepted'] >= math.ceil(.8 * bootstrap) and boot['pc_interval'] is not None and
        rates[0] < boot['pc_interval'][0] < boot['pc_interval'][1] < rates[-1] and
        H6.NU_BOUNDS[0] < boot['nu_interval'][0] < boot['nu_interval'][1] < H6.NU_BOUNDS[1])
    variants = [('drop_smallest', [r for r in selected if r['steps'] != sizes[0]], 2, False),
                ('narrow_window', [r for r in selected if rates[0] < r['p'] < rates[-1]], 2, False),
                ('cubic', selected, 3, False), ('inverse_L_correction', selected, 2, True)]
    sensitivity = []
    for name, subset, degree, correction in variants:
        enough = len({r['steps'] for r in subset}) >= 4 and len({r['p'] for r in subset}) >= 5
        fit = H6.candidate_fit(subset, degree, correction) if enough else None
        pc_tol = max(.1 * (rates[-1] - rates[0]), np.ptp(boot['pc_interval']) if boot['pc_interval'] else 0.)
        nu_tol = max(.25 * candidate['nu'], np.ptp(boot['nu_interval']) if boot['nu_interval'] else 0.) if candidate else None
        stable = bool(candidate and fit and abs(fit['pc'] - candidate['pc']) <= pc_tol and abs(fit['nu'] - candidate['nu']) <= nu_tol)
        passed = bool(enough and stable and all(H6.numerical_guards(fit).values()))
        guards['sensitivity_' + name] = passed
        sensitivity.append(dict(name=name, candidate=fit, enough=enough, stable=stable, passed=passed,
                                pc_tolerance=float(pc_tol), nu_tolerance=float(nu_tol) if nu_tol else None))
    accepted = all(guards.values())
    return dict(status='accepted' if accepted else 'diagnostic_only',
                pc=candidate['pc'] if accepted else None, nu=candidate['nu'] if accepted else None,
                candidate=candidate, sizes=sizes, rates=rates, guards=guards,
                failed_guards=[k for k, v in guards.items() if not v],
                crossings=crossings, crossing_common_brackets=clusters, bootstrap=boot,
                sensitivity=sensitivity, guard_policy=H6.GUARD_POLICY)


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_figures(points, marginals, differences, fit, output):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/toponoise-h6-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from plotting_theme import paper_theme, save_figure, style_axes
    sizes = sorted({r['steps'] for r in points})
    # Existing study colors stay fixed; other allowed sizes have a stable color
    # determined by L itself, independent of the sizes present in a particular plot.
    colors = {L: COLORS.get(L, plt.get_cmap('turbo')(.15 + .7 * ((L * .6180339887498949) % 1))) for L in sizes}
    rates = sorted({r['p'] for r in points})
    def curves(ax, arm):
        for L in sizes:
            rows = [r for r in marginals if r['arm'] == arm and r['L'] == L]
            y = np.array([r['rate'] for r in rows])
            ax.errorbar([r['p'] for r in rows], y, yerr=np.maximum(0, [y - np.array([r['ci_low'] for r in rows]),
                np.array([r['ci_high'] for r in rows]) - y]), marker='o', color=colors[L], label=fr'$L={L}$', elinewidth=.85)
        style_axes(ax)
        ax.set(xlabel=r'Physical Z fault probability $p$', ylabel='Logical Z state failure rate', ylim=(0, 1))
        ax.margins(x=.04)
    with paper_theme():
        fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.8), sharex=True, sharey=True)
        for ax, arm, title in zip(axes.flat, ARMS, LABELS):
            curves(ax, arm)
            ax.set_title(title)
            ax.legend(loc='upper left', ncol=2, fontsize=9)
        fig.subplots_adjust(left=.08, right=.98, bottom=.10, top=.94, wspace=.23, hspace=.28)
        fig.text(.08, .025, 'H+6 · Z-only · T = L · same fault records for all four arms · 95% Wilson intervals', fontsize=9)
        save_figure(fig, output / 'z-one-step-four-protocols')
        plt.close(fig)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), sharex=True, sharey=True)
        for ax, first, title in zip(axes, ('final_only', 'seed_only', 'window_complete'),
                                   ('One-step − final-only', 'One-step − seed-only', 'One-step − complete-window')):
            for L in sizes:
                rows = [r for r in differences if r['first'] == first and r['L'] == L]
                y = np.array([r['difference'] for r in rows])
                ax.errorbar([r['p'] for r in rows], y, yerr=np.maximum(0, [y - np.array([r['ci_low'] for r in rows]),
                    np.array([r['ci_high'] for r in rows]) - y]), marker='o', color=colors[L], label=fr'$L={L}$', elinewidth=.85)
            ax.axhline(0, color='0.4', ls='--', lw=.9)
            style_axes(ax)
            ax.set(title=title, xlabel=r'Physical Z fault probability $p$')
        axes[0].set_ylabel('Difference in failure probability')
        axes[0].legend(loc='best', ncol=2, fontsize=8)
        fig.subplots_adjust(left=.065, right=.98, bottom=.19, top=.89, wspace=.20)
        fig.text(.065, .035, 'Negative favors one-step · approximate pointwise 95% paired intervals (Bonferroni–Wilson); no multiple-comparison adjustment', fontsize=9)
        save_figure(fig, output / 'z-one-step-paired-differences')
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8.2, 5.5))
        curves(ax, 'one_step')
        ax.set_title('H+6 · Causal one-step Z feedback')
        ax.legend(loc='upper right', ncol=2, fontsize=9)
        candidate = fit['candidate'] if fit else None
        if candidate and not candidate['bound_hits']:
            # Show the rejected candidate transparently; plotting it never
            # changes fit acceptance or publishes a threshold estimate.
            inset = ax.inset_axes([.12, .69, .38, .23])
            selected = [r for r in marginals if r['arm'] == 'one_step' and r['p'] in fit['rates']]
            xs = []
            for L in sizes:
                rows = [r for r in selected if r['L'] == L]
                x = [(r['p'] - candidate['pc']) * L**(1/candidate['nu']) for r in rows]
                xs.extend(x)
                y = np.array([r['rate'] for r in rows])
                inset.errorbar(x, y, yerr=np.maximum(0, [y - np.array([r['ci_low'] for r in rows]),
                    np.array([r['ci_high'] for r in rows]) - y]), marker='o', color=colors[L],
                    ms=2.5, lw=.9, elinewidth=.55, capsize=1)
            grid = np.linspace(min(xs), max(xs), 300)
            xnorm = grid / (candidate['p_scale'] * candidate['L_reference']**(1/candidate['nu']))
            inset.plot(grid, expit(np.polynomial.polynomial.polyval(xnorm, candidate['coefficients'])),
                       color='0.2', lw=1, ls='-' if fit['status'] == 'accepted' else '--')
            inset.set(xlabel=r'$(p-p_c)L^{1/\nu}$', ylabel='Logical error')
            inset.xaxis.set_major_locator(MaxNLocator(4))
            inset.set_title('Accepted fit' if fit['status'] == 'accepted' else 'Rejected fit: diagnostic only', fontsize=8)
            style_axes(inset, inset=True)
            if fit['status'] == 'accepted':
                ax.axvline(candidate['pc'], color='0.4', ls='--', lw=.9, zorder=0)
        text = ('Threshold not established' if not fit or fit['status'] != 'accepted' else
                fr'$p_c={fit["pc"]:.5f}$; $\nu={fit["nu"]:.2f}$')
        ax.text(.97, .16, text, transform=ax.transAxes, ha='right', fontsize=9)
        fig.subplots_adjust(left=.12, right=.97, bottom=.16, top=.90)
        fig.text(.12, .025, '95% Wilson bars · complete seed correction + retiring-sheet transport · T = L', fontsize=9)
        save_figure(fig, output / 'z-one-step-local-fit')
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', action='append', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--selection', type=Path)
    parser.add_argument('--bootstrap', type=int, default=200)
    parser.add_argument('--seed', type=int, default=2026101514)
    args = parser.parse_args()
    points, inputs = load_points(args.input)
    selection = json.loads(args.selection.read_text()) if args.selection else None
    selection_validation = (validate_selection(points, inputs, selection, args.selection,
                            bootstrap=args.bootstrap, seed=args.seed) if selection is not None else None)
    marginals, differences = [], []
    for row in points:
        for arm in ARMS:
            k, n = row[arm + '_failures'], row['shots']
            lo, hi = H6.wilson_interval(k, n)
            marginals.append(dict(L=row['steps'], p=row['p'], arm=arm, failures=k, shots=n,
                                  rate=k/n, ci_low=lo, ci_high=hi))
        for first in ('final_only', 'seed_only', 'window_complete'):
            differences.append(dict(L=row['steps'], p=row['p'], shots=row['shots'],
                                    **paired(row['joint'], first, 'one_step')))
    fit = scaling_diagnostics(points, selection['fit_rates'], args.bootstrap, args.seed) if selection and selection.get('fit_rates') else None
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / 'marginals.csv', marginals)
    write_csv(output / 'paired-differences.csv', differences)
    report = dict(arms=list(ARMS), inputs=inputs, selection=selection, selection_validation=selection_validation,
                  total_fault_records=sum(r['shots'] for r in points),
                  total_protocol_outcomes=4 * sum(r['shots'] for r in points),
                  points=points, marginals=marginals, paired_differences=differences, fit=fit,
                  analysis_source_sha256=H6.sha256(__file__),
                  interpretation='Paired fixed-p comparisons are distinct from threshold inference. Seed-only isolates intermediate feedback. No Z/X threshold comparison is implied.')
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    plot_figures(points, marginals, differences, fit, output)
    print(json.dumps(dict(points=len(points), records=report['total_fault_records'],
                          fit_status=fit['status'] if fit else 'not_requested'), indent=2))


if __name__ == '__main__':
    main()
