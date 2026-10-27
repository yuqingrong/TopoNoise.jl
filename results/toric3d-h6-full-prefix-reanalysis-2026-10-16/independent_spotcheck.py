#!/usr/bin/env python3
"""Read-only, independent arithmetic checks of the saved exploratory analyses.

Does not import either analysis implementation or rerun the full bootstrap or
profiles. Three fixed-(pc, nu), convex nuisance fits provide bounded spot checks.
"""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.special import expit, gammaln, xlogy
from scipy.optimize import minimize
from scipy.stats import binomtest, chi2

HERE = Path(__file__).resolve().parent
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def close(actual, expected, message, tolerance=1e-8):
    require(np.allclose(actual, expected, atol=tolerance, rtol=0), message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def merge(regions):
    result = []
    for left, right in sorted(regions):
        if result and left <= result[-1][1]:
            result[-1][1] = max(result[-1][1], right)
        else:
            result.append([left, right])
    return result


def shared_regions(rates, intervals, pairs):
    """Solve all affine overlap inequalities independently on each grid segment."""
    result = []
    for index, (left, right) in enumerate(zip(rates, rates[1:])):
        candidates = [0., 1.]
        constraints = []
        for a, b in pairs:
            for small, large in ((a, b), (b, a)):
                f0 = intervals[small, rates[index]][0] - intervals[large, rates[index]][1]
                f1 = intervals[small, rates[index+1]][0] - intervals[large, rates[index+1]][1]
                constraints.append((f0, f1))
                if f0 != f1:
                    root = -f0/(f1-f0)
                    if 0 <= root <= 1:
                        candidates.append(root)
        candidates = sorted(set(candidates))
        valid = lambda t: all(a+(b-a)*t <= 1e-14 for a, b in constraints)
        for a, b in zip(candidates, candidates[1:]):
            if valid((a+b)/2):
                result.append([left+(right-left)*a, left+(right-left)*b])
        for t in candidates:
            if valid(t):
                result.append([left+(right-left)*t, left+(right-left)*t])
    return merge(result)


def deviance(counts, shots, probability):
    obs = counts/shots
    return float(2*np.sum(xlogy(counts, obs/probability) +
                         xlogy(shots-counts, (1-obs)/(1-probability))))


def old_guard(counts, points, pc):
    sizes = sorted({row['steps'] for row in points})
    rates = sorted({row['p'] for row in points})
    table = {(row['steps'], row['p']): count/row['shots'] for row, count in zip(points, counts)}
    common = [[rates[0], rates[-1]]]
    for smaller, larger in zip(sizes, sizes[1:]):
        diff = [table[larger, p]-table[smaller, p] for p in rates]
        brackets = [[rates[i], rates[i+1]] for i in range(len(rates)-1) if diff[i]*diff[i+1] < 0]
        brackets += [[rates[i-1], rates[i+1]] for i in range(1, len(rates)-1)
                     if diff[i] == 0 and diff[i-1]*diff[i+1] < 0]
        common = merge([[max(a, c), min(b, d)] for a, b in common for c, d in brackets
                        if max(a, c) <= min(b, d)])
    return not (len(common) == 1 and common[0][0] <= pc <= common[0][1])


def main():
    snapshot = json.loads((HERE/'input.json').read_text())
    freeze = json.loads((HERE/'plan-freeze.json').read_text())
    require(sha(HERE/'input.json') == freeze['input_sha256'], 'frozen input hash')
    require(sha(HERE/'analysis-plan.md') == freeze['plan_sha256'], 'frozen plan hash')
    source = Path(snapshot['source_directory'])
    require(sha(source/'analysis/report.json') == snapshot['source_report_sha256'], 'original report hash')
    require(sha(source/'ARTIFACTS.sha256.json') == snapshot['source_manifest_sha256'], 'original manifest hash')
    original = json.loads((source/'analysis/report.json').read_text())
    require(len(snapshot['points']) == len(original['points']) == 70, '70 fixed cells')
    for row, old in zip(snapshot['points'], original['points']):
        require({key: row[key] for key in old} == old, 'unchanged source point')
        require(row['failures'] == row['none_failures'] == row['full_prefix_failures'], 'failure aliases')
    require(sum(row['shots'] for row in snapshot['points']) == 70000, '70,000 existing records')

    x = json.loads((HERE/'x-report.json').read_text())
    z = json.loads((HERE/'z-report.json').read_text())
    require(z.get('completed_at_unix') is not None, 'Z report completed')
    require(x['provenance']['analysis_sha256'] == sha(HERE/'x_analysis.py'), 'X source binding')
    require(z['source_sha256'] == sha(HERE/'z_models.py'), 'Z source binding')
    require(x['provenance']['input_sha256'] == z['input_sha256'] == freeze['input_sha256'], 'same input')
    require(x['provenance']['plan_sha256'] == z['plan_sha256'] == freeze['plan_sha256'], 'same plan')
    require(freeze['fixed_before_followup_fits_at_unix'] < x['provenance']['started_at_unix'], 'X freeze chronology')
    require(freeze['fixed_before_followup_fits_at_unix'] < z['created_at_unix'], 'Z freeze chronology')

    points = [row for row in snapshot['points'] if row['channel'] == 'x_only']
    sizes = sorted({row['steps'] for row in points})
    rates = sorted({row['p'] for row in points})
    lookup = {(row['steps'], row['p']): row for row in points}
    intervals = {}
    for saved in x['crossing']['cells']:
        original_cell = lookup[saved['L'], saved['p']]
        require((saved['failures'], saved['shots']) ==
                (original_cell['failures'], original_cell['shots']), 'X cell counts from frozen input')
        interval = binomtest(saved['failures'], saved['shots']).proportion_ci(1-.05/35, method='exact')
        close([saved['ci_low'], saved['ci_high']], [interval.low, interval.high], 'independent CP interface', 5e-12)
        intervals[saved['L'], saved['p']] = (interval.low, interval.high)
    adjacent = shared_regions(rates, intervals, list(zip(sizes, sizes[1:])))
    all_pairs = shared_regions(rates, intervals, list(itertools.combinations(sizes, 2)))
    close(x['crossing']['primary_adjacent']['components'], adjacent, 'primary interval geometry', 5e-12)
    close(x['crossing']['secondary_all_pairs']['components'], all_pairs, 'all-pair interval geometry', 5e-12)
    n = np.array([row['shots'] for row in points])
    k = np.array([row['failures'] for row in points])
    common = x['common_model']; independent = x['independent_size_model']
    xcoord = np.array([(row['p']-common['pc'])/common['p_scale'] *
                      (row['steps']/common['L_reference'])**(1/common['nu']) for row in points])
    probability = expit(sum(coefficient*xcoord**power for power, coefficient in enumerate(common['coefficients'])))
    close(probability, common['fitted_probabilities'], 'X common predictions')
    close(deviance(k, n, probability), common['deviance'], 'X common deviance')
    by_size = {row['L']: row['coefficients'] for row in independent['curves']}
    t = np.array([2*(row['p']-rates[0])/(rates[-1]-rates[0])-1 for row in points])
    alt_probability = expit([sum(c*t[i]**power for power, c in enumerate(by_size[row['steps']]))
                             for i, row in enumerate(points)])
    close(alt_probability, independent['fitted_probabilities'], 'independent curve predictions')
    close(deviance(k, n, alt_probability), independent['deviance'], 'independent curve deviance')
    close(x['observed_likelihood_ratio'], common['deviance']-independent['deviance'], 'observed LR')
    for size in sizes:
        mask = np.array([row['steps'] == size for row in points])
        score = np.column_stack([np.ones(mask.sum()), t[mask], t[mask]**2]).T @ (n[mask]*alt_probability[mask]-k[mask])
        require(np.max(np.abs(score))/n[mask].sum() < 1e-9, 'independent convex likelihood score')

    bootstrap = x['bootstrap']; records = bootstrap['records']
    require(len(records) == bootstrap['replicates'] == 1000, 'all 1000 bootstrap records')
    require([row['replicate'] for row in records] == list(range(1000)), 'bootstrap sequence')
    draws = np.random.default_rng(bootstrap['seed']).binomial(n, probability, size=(1000, len(points)))
    for row in records:
        require(row['simulated_failures'] == draws[row['replicate']].tolist(), 'saved bootstrap counts match regenerated stream')
    valid = [row for row in records if row['common_valid']]
    joint = [row for row in valid if row['alternative_valid']]
    for row in valid:
        require(row['old_crossing_rejected'] == old_guard(draws[row['replicate']], points, row['pc']), 'old crossing guard from regenerated counts')
    for label, subset, field, observed in (
        ('goodness_of_fit', valid, 'common_deviance', common['deviance']),
        ('likelihood_ratio', joint, 'likelihood_ratio', x['observed_likelihood_ratio'])):
        summary = bootstrap[label]
        tail = sum(row[field] >= observed-1e-10 for row in subset)
        require(summary['tail_count'] == tail and summary['valid_replicates'] == len(subset), label+' denominator')
        close(summary['p_plus_one'], (tail+1)/(len(subset)+1), label+' plus one')
        interval = binomtest(tail, len(subset)).proportion_ci(.95, method='exact')
        close(summary['mc_binomial_interval'], [interval.low, interval.high], label+' MC CI', 5e-12)
    bootstrap_spots = []
    # A bounded independent fit of only three linear coefficients at three saved
    # (pc, nu) pairs verifies that regenerated counts underlie the reported NLLs.
    # This does not rerun the nonlinear common-model bootstrap or any profiles.
    for index in (0, 499, 999):
        row = records[index]; counts = draws[index]
        coord = np.array([(p['p']-row['pc'])/common['p_scale'] *
                          (p['steps']/common['L_reference'])**(1/row['nu']) for p in points])
        design = np.column_stack([np.ones(len(coord)), coord, coord**2])
        def objective(coefficients):
            eta = design@coefficients
            return (float(np.sum(n*np.logaddexp(0., eta)-counts*eta)/n.sum()),
                    design.T@(n*expit(eta)-counts)/n.sum())
        result = minimize(objective, common['coefficients'], jac=True, method='BFGS',
                          options={'gtol': 1e-10, 'maxiter': 500})
        require(np.max(np.abs(objective(result.x)[1])) < 1e-8, 'bounded convex bootstrap spot fit')
        computed = deviance(counts, n, expit(design@result.x))
        close(computed, row['common_deviance'], 'bootstrap saved deviance spot check', 1e-5)
        bootstrap_spots.append({'replicate': index, 'reported_deviance': row['common_deviance'],
                                'independent_deviance': computed})
    old_count = sum(row['old_crossing_rejected'] for row in valid)
    require(bootstrap['old_raw_crossing_rule']['tail_count'] == old_count, 'old guard count')
    close(bootstrap['old_raw_crossing_rule']['rejection_frequency'], old_count/len(valid), 'old guard rate')
    for field in ('pc', 'nu'):
        close(bootstrap['conditional_percentile_intervals'][field], np.quantile([row[field] for row in valid], [.025, .975]), field+' percentile interval')
    require(bootstrap['bound_hit_replicates'] == sum(bool(row['bound_hits']) for row in valid), 'bound flags retained')
    require(bootstrap['nonmonotone_replicates'] == sum(not row['monotone'] for row in valid), 'monotonicity flags retained')

    require(len(z['fits']) == 36 and len(z['profiles']) == 4, 'fixed Z model/profile counts')
    models = {'quadratic', 'cubic', 'corr_w0p5', 'corr_w1', 'corr_w2', 'corr_linear_w1'}
    require({(fit['min_size'], fit['window'], fit['model']) for fit in z['fits']} ==
            set(itertools.product((3, 5, 7), ('full', 'narrow'), models)), 'complete 36-model matrix')
    fits = {fit['fit_id']: fit for fit in z['fits']}
    z_lookup = {(row['steps'], row['p']): row for row in snapshot['points'] if row['channel'] == 'z_only'}
    for fit in fits.values():
        rows = fit['predictions']; n = np.array([row['shots'] for row in rows]); k = np.array([row['failures'] for row in rows])
        expected_keys = {(size, p) for size, p in z_lookup if size >= fit['min_size'] and fit['p_min'] <= p <= fit['p_max']}
        require({(row['size'], row['p']) for row in rows} == expected_keys, 'Z prescribed subset cells')
        for row in rows:
            original_cell = z_lookup[row['size'], row['p']]
            require((row['failures'], row['shots']) ==
                    (original_cell['full_prefix_failures'], original_cell['shots']), 'Z counts from frozen input')
        size = np.array([row['size'] for row in rows]); p = np.array([row['p'] for row in rows])
        coordinate = (p-fit['pc'])/fit['p_scale']*(size/fit['L_reference'])**(1/fit['nu'])
        eta = sum(c*coordinate**power for power, c in enumerate(fit['coefficients']))
        if fit['correction_omega'] is not None:
            c = fit['correction_coefficients']
            eta += (c[0] + (c[1]*coordinate if len(c) == 2 else 0))*(fit['L_reference']/size)**fit['correction_omega']
        probability = expit(eta)
        nll = float(np.sum(n*np.logaddexp(0., eta)-k*eta))
        ll = float(np.sum(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)))-nll
        count = fit['n_parameters']; cells = len(rows)
        close(probability, fit['fitted_probabilities'], 'Z fitted probabilities')
        close(fit['nll'], nll, 'Z total NLL')
        close(fit['log_likelihood'], ll, 'Z full binomial log likelihood')
        close(fit['deviance'], deviance(k, n, probability), 'Z deviance')
        require(fit['degrees_of_freedom'] == cells-count, 'Z df')
        close(fit['deviance_p_value'], chi2.sf(fit['deviance'], cells-count), 'Z asymptotic goodness')
        close(fit['aicc'], -2*ll+2*count+2*count*(count+1)/(cells-count-1), 'Z AICc')
        floor = min(other['aicc'] for other in fits.values() if (other['min_size'], other['window']) == (fit['min_size'], fit['window']))
        close(fit['delta_aicc_same_subset'], fit['aicc']-floor, 'Z within-subset AICc delta')
    for profile in z['profiles']:
        fit = fits[profile['fit_id']]
        predictions = fit['predictions']
        sizes = np.array([row['size'] for row in predictions]); rates = np.array([row['p'] for row in predictions])
        shots = np.array([row['shots'] for row in predictions]); counts = np.array([row['failures'] for row in predictions])
        require(profile['minimum_agreement'], 'profile/free minimum agreement')
        close(profile['lr_cutoff'], chi2.ppf(.95, 1), 'profile cutoff')
        successful = [row for row in profile['evaluations'] if row['optimizer_success']]
        require(len(profile['evaluations']) >= 301, 'profile grid resolution')
        for row in successful:
            close(row['lr'], 2*(row['nll']-fit['nll']), 'profile total LR', 2e-8)
            require(row['lr'] >= -4e-4, 'profile cannot improve free fit materially')
            parameters = row['parameters']; nu = math.exp(parameters[1])
            close(parameters[0], (row['pc']-fit['p_min'])/fit['p_scale'], 'profile fixed pc')
            coord = (rates-row['pc'])/fit['p_scale']*(sizes/fit['L_reference'])**(1/nu)
            eta = parameters[2]+parameters[3]*coord+parameters[4]*coord**2
            if fit['model'] == 'corr_w1':
                eta += parameters[5]*fit['L_reference']/sizes
            independent_nll = float(np.sum(shots*np.logaddexp(0., eta)-counts*eta))
            close(row['nll'], independent_nll, 'profile independent likelihood reconstruction', 2e-8)
        for region in profile['regions']:
            for side in ('lower', 'upper'):
                if region[side+'_bound_truncated'] or region[side+'_numerical_gap']:
                    continue
                endpoint = min(successful, key=lambda row: abs(row['pc']-region[side]))
                close(endpoint['pc'], region[side], 'saved refined endpoint', 2e-12)
                close(endpoint['lr'], profile['lr_cutoff'], 'refined endpoint cutoff', 2e-5)

    print(json.dumps({'status': 'passed', 'assertions': checks,
        'scope': 'Independent saved-output arithmetic, fixed-input/source bindings, interval geometry, 1000 regenerated bootstrap old guards, three fixed-(pc,nu) convex nuisance fits; no analysis imports or full simulation/bootstrap/profile reruns.',
        'input_sha256': sha(HERE/'input.json'), 'x_report_sha256': sha(HERE/'x-report.json'),
        'z_report_sha256': sha(HERE/'z-report.json'), 'review_script_sha256': sha(Path(__file__)),
        'runtime': {'python_executable': sys.executable, 'python_version': sys.version,
                    'numpy': np.__version__, 'scipy': scipy.__version__,
                    'thread_environment': {key: os.environ.get(key) for key in
                        ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}},
        'x_primary_components': adjacent, 'x_all_pair_components': all_pairs,
        'x_bootstrap_valid_common': len(valid), 'x_bootstrap_valid_both': len(joint),
        'bounded_convex_bootstrap_spot_fits': bootstrap_spots,
        'z_profile_numerical_failures': {row['fit_id']: len(row['numerical_failures']) for row in z['profiles']}}, indent=2))


if __name__ == '__main__':
    main()
