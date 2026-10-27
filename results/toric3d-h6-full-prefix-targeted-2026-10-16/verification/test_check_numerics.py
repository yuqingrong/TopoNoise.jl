"""Small mathematical fixtures for the independent saved-output audit."""
import importlib.util
import copy
import math
from pathlib import Path
import unittest

import numpy as np


class IndependentNumericsTests(unittest.TestCase):
    def checker(self):
        path = Path(__file__).with_name("check_numerics.py")
        self.assertTrue(path.is_file(), "independent numerical checker is missing")
        spec = importlib.util.spec_from_file_location("independent_numerics", path)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        return module

    def test_unequal_binomial_counts_and_saturated_deviance(self):
        a = self.checker()
        value = a.likelihood_stats(np.array([5., 10.]), np.array([10., 20.]), np.zeros(2))
        self.assertAlmostEqual(value["nll"], 30 * math.log(2))
        self.assertAlmostEqual(value["log_likelihood"], math.log(252 * 184756) - 30 * math.log(2))
        self.assertAlmostEqual(value["deviance"], 0.)
        with self.assertRaises(a.AuditError):
            a.close(value["nll"] + 1, 30 * math.log(2), "corrupt likelihood")

    def test_linear_correction_uses_size_and_scaled_coordinate(self):
        a = self.checker()
        eta = a.evaluate_parameters([.5, 0., .1, .2, .3, .4, .5], np.array([1., 2.]),
                                    np.array([.05, .1]), (0., .1), 1., "corr_linear_w1")
        np.testing.assert_allclose(eta, [.5, 1.05], atol=1e-14)

    def test_boundary_and_monotonicity_flags_are_recomputed(self):
        a = self.checker()
        result = a.parameter_diagnostics([.001, 0., 0., -1., 0.], "quadratic", [5, 11], (0., .1), 5.)
        self.assertEqual(result["bound_hits"], ["pc"])
        self.assertFalse(result["monotone_increasing"])
        with self.assertRaises(a.AuditError):
            a.parameter_diagnostics([1.1, 0., 0., 1., 0.], "quadratic", [5, 11], (0., .1), 5.)

    def test_refined_profile_endpoint_must_meet_likelihood_cutoff(self):
        a = self.checker()
        cutoff = 3.841458820694124
        profile = {"lr_cutoff": cutoff, "pc_bounds": [.02, .05], "regions": [
            {"lower": .03, "upper": .04, "lower_bound_truncated": False, "upper_bound_truncated": False,
             "lower_numerical_gap": False, "upper_numerical_gap": False, "refinement_failures": []}]}
        evaluations = [{"pc": .03, "lr": cutoff}, {"pc": .035, "lr": 0.}, {"pc": .04, "lr": cutoff}]
        a.validate_regions(profile, evaluations)
        evaluations[-1]["lr"] = 1.
        with self.assertRaises(a.AuditError):
            a.validate_regions(profile, evaluations)

    def test_complete_fit_record_checks_likelihood_and_source_counts(self):
        a = self.checker()
        rows = [{"steps":L,"p":p,"shots":n,"failures":n//2} for (L,p),n in
                zip([(5,.025),(5,.0325),(5,.04),(11,.025),(11,.0325),(11,.04)], [10,20,30,40,50,60])]
        nll = 210*math.log(2)
        ll = sum(math.log(math.comb(r["shots"],r["failures"])) for r in rows)-nll
        fit = {"fit_id":"fixture","model":"quadratic","parameters":[.5,0.,0.,0.,0.],
               "L_reference":math.sqrt(55),"p_min":.025,"p_max":.04,"p_scale":.015,"pc":.0325,"nu":1.,
               "parameter_names":["pc","nu","a0","a1","a2"],
               "parameter_bounds":[[0.,1.],[math.log(.3),math.log(5.)],[-15.,15.],[-30.,30.],[-30.,30.]],
               "coefficients":[0.,0.,0.],"correction_omega":None,"correction_coefficients":[],
               "bound_hits":[],"monotone_increasing":False,"derivative_min":0.,"nll":nll,"log_likelihood":ll,
               "deviance":0.,"fitted_probabilities":[.5]*6,"standardized_residuals":[0.]*6,
               "largest_absolute_residual":0.,"n_cells":6,"n_parameters":5,"degrees_of_freedom":1,
               "diagnostic_three_sizes_only":False,"deviance_p_value":1.,"aic":-2*ll+10,"aicc":None,"aicc_sample_size":6,
               "predictions":[{"size":r["steps"],"p":r["p"],"shots":r["shots"],"failures":r["failures"],
                               "predicted":.5,"observed":.5,"pearson_residual":0.,"deviance_residual":0.} for r in rows]}
        a.audit_fit(fit,rows,[.025,.04])
        wrong = copy.deepcopy(fit); wrong["nll"] += 1
        with self.assertRaises(a.AuditError): a.audit_fit(wrong,rows,[.025,.04])
        wrong = copy.deepcopy(fit); wrong["predictions"][0]["shots"] = 12
        with self.assertRaises(a.AuditError): a.audit_fit(wrong,rows,[.025,.04])


if __name__ == "__main__":
    unittest.main()
