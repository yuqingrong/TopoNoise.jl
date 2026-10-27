import importlib.util
import json
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
try:
    import z_models as z
except ModuleNotFoundError:
    z = None


class Availability(unittest.TestCase):
    def test_z_model_analysis_api_exists(self):
        self.assertIsNotNone(z, 'The fixed-grid Z likelihood analysis is not implemented')


@unittest.skipIf(z is None, 'analysis API is absent')
class Models(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.input = json.loads((HERE/'input.json').read_text())
        cls.data = z.select_data(cls.input['points'], 3, 'full')

    def test_analytic_gradient_for_all_six_models(self):
        for model in z.MODELS:
            with self.subTest(model=model):
                theta = np.array([.57, math.log(1.17), -.3, 2.1, -.65] +
                    ([.17] if model == 'cubic' else [.21, -.13] if model == 'corr_linear_w1'
                     else [.21] if model.startswith('corr_') else []))
                loss, gradient = z.objective(theta, self.data, model)
                numerical = np.empty_like(theta)
                for j in range(len(theta)):
                    plus, minus = theta.copy(), theta.copy()
                    plus[j] += 1e-6; minus[j] -= 1e-6
                    numerical[j] = (z.objective(plus, self.data, model)[0] -
                                    z.objective(minus, self.data, model)[0]) / 2e-6
                self.assertTrue(np.isfinite(loss))
                np.testing.assert_allclose(gradient, numerical, rtol=2e-5, atol=2e-9)

    def test_correction_models_contain_leading_quadratic(self):
        base = np.array([.6, math.log(1.2), -.2, 2., -.7])
        expected, _ = z.eta_jacobian(base, self.data, 'quadratic')
        for model in z.MODELS:
            theta = np.r_[base, np.zeros(z.parameter_count(model)-5)]
            eta, jacobian = z.eta_jacobian(theta, self.data, model)
            np.testing.assert_allclose(eta, expected, rtol=0, atol=1e-14)
            self.assertEqual(jacobian.shape, (len(expected), len(theta)))

    def test_original_quadratic_and_omega1_fits_are_reproduced(self):
        original = self.input['original_fits']['z_only']
        targets = {'quadratic': original['candidate'], 'corr_w1': next(
            item['candidate'] for item in original['sensitivity']
            if item['name'] == 'inverse_L_correction')}
        for model, target in targets.items():
            result = z.fit_model(self.data, model)
            self.assertTrue(result['optimizer_success'])
            self.assertAlmostEqual(result['pc'], target['pc'], delta=2e-7)
            self.assertAlmostEqual(result['nu'], target['nu'], delta=5e-5)
            self.assertAlmostEqual(result['deviance'], target['deviance'], delta=2e-5)
            np.testing.assert_allclose(result['fitted_probabilities'],
                                       target['fitted_probabilities'], atol=3e-6, rtol=0)

    def test_known_synthetic_quadratic_and_correction_parameters(self):
        for model in ('quadratic', 'corr_w1'):
            theta = np.array([.53, math.log(1.1), -.5, 2., -.4] +
                             ([.25] if model == 'corr_w1' else []))
            eta, _ = z.eta_jacobian(theta, self.data, model)
            synthetic = dict(self.data)
            synthetic['shots'] = np.full(len(eta), 100000000., dtype=float)
            synthetic['failures'] = np.rint(synthetic['shots'] * z.expit(eta))
            fitted = z.fit_model(synthetic, model)
            self.assertTrue(fitted['optimizer_success'])
            self.assertAlmostEqual(fitted['pc'], .025+.015*.53, delta=2e-6)
            self.assertAlmostEqual(fitted['nu'], 1.1, delta=5e-4)
            self.assertLess(fitted['deviance'], .001)

    def test_monotonicity_checks_correction_slope_per_size(self):
        theta = np.array([.5, 0., 0., 1., 0., 0., -1.])
        diagnostic = z.monotonicity(theta, self.data, 'corr_linear_w1')
        self.assertFalse(diagnostic['monotone_increasing'])
        self.assertLess(diagnostic['derivative_min'], 0)

    def test_profile_reoptimizes_nuisance_and_agrees_at_free_minimum(self):
        fit = z.fit_model(self.data, 'quadratic')
        point = z.fit_fixed_pc(self.data, 'quadratic', fit['pc'], [fit['parameters']])
        self.assertTrue(point['optimizer_success'])
        self.assertAlmostEqual(point['nll'], fit['nll'], delta=2e-5)
        off = z.fit_fixed_pc(self.data, 'quadratic', fit['pc']+.0005, [fit['parameters']])
        self.assertGreater(off['nll'], fit['nll'])
        self.assertGreater(np.linalg.norm(np.array(off['parameters'])[1:] -
                                          np.array(fit['parameters'])[1:]), 1e-5)

    def test_failed_polish_cannot_replace_an_equivalent_converged_solution(self):
        success = SimpleNamespace(fun=1., x=np.array([.5]), jac=np.array([0.]), success=True, message='converged')
        failed = SimpleNamespace(fun=1.-1e-15, x=np.array([.5]), jac=np.array([1e-10]), success=False, message='ABNORMAL')
        with patch.object(z, 'minimize', side_effect=[success, failed]):
            result, diagnostic = z._optimize(lambda x:(1.,np.zeros(1)), [(0.,1.)], [np.array([.5])])
        self.assertTrue(result.success)
        self.assertFalse(diagnostic['polish_success'])

    def test_profile_rejects_a_minimum_below_the_free_fit(self):
        fit = z.fit_model(self.data, 'quadratic')
        fit['fit_id'] = 'test-quadratic'
        fit['nll'] += .01
        with self.assertRaisesRegex(RuntimeError, 'improves the free fit'):
            z.profile_pc(self.data, 'quadratic', fit, grid_count=5)

    def test_connected_regions_include_all_islands_and_bound_truncation(self):
        curve = lambda p: min(((p-.2)/.1)**2, ((p-.75)/.1)**2)
        grid = np.linspace(0, 1, 101)
        regions = z.connected_regions(grid, [curve(p) for p in grid], 1., curve)
        self.assertEqual(len(regions), 2)
        np.testing.assert_allclose([[r['lower'],r['upper']] for r in regions],
                                   [[.1,.3],[.65,.85]], atol=1e-8)
        bounded = z.connected_regions(grid, grid**2, .04, lambda p:p*p)
        self.assertTrue(bounded[0]['lower_bound_truncated'])
        self.assertAlmostEqual(bounded[0]['upper'], .2, delta=1e-8)
        gaps = z.connected_regions(np.array([0., .25, .5, .75, 1.]),
                                   [0.,0.,float('nan'),0.,0.], 1., lambda p:0.)
        self.assertEqual(len(gaps), 2)
        self.assertTrue(gaps[0]['upper_numerical_gap'])
        self.assertTrue(gaps[1]['lower_numerical_gap'])

    def test_failed_root_refinement_is_retained_as_a_numerical_gap(self):
        curve = lambda p: float('nan') if .2<p<.3 else ((p-.5)/.25)**2
        regions = z.connected_regions(np.array([0.,.5,1.]), [4.,0.,4.], 1., curve)
        self.assertEqual(len(regions),1)
        self.assertTrue(regions[0]['lower_numerical_gap'])
        self.assertEqual(len(regions[0]['refinement_failures']),1)
        self.assertAlmostEqual(regions[0]['upper'],.75,delta=1e-8)

    def test_aicc_uses_binomial_cells_and_full_binomial_likelihood(self):
        fitted = z.fit_model(self.data, 'quadratic')
        n, k = 35, 5
        self.assertEqual(fitted['aicc_sample_size'], n)
        expected = -2*fitted['log_likelihood'] + 2*k + 2*k*(k+1)/(n-k-1)
        self.assertAlmostEqual(fitted['aicc'], expected, places=8)
        self.assertIn('heuristic', fitted['aicc_interpretation'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
