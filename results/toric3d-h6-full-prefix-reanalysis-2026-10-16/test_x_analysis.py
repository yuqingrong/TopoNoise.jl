import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.special import expit
from scipy.stats import binom

PATH=Path(__file__).with_name('x_analysis.py')
M=None
if PATH.exists():
    spec=importlib.util.spec_from_file_location('x_reanalysis',PATH)
    M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)


class XAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(M,'X uncertainty-aware analysis is not implemented')

    def test_exact_binomial_endpoints_and_known_case(self):
        lo,hi=M.clopper_pearson(0,10,.05)
        self.assertEqual(lo,0.)
        self.assertAlmostEqual(hi,1-.025**.1)
        lo,hi=M.clopper_pearson(10,10,.05)
        self.assertAlmostEqual(lo,.025**.1)
        self.assertEqual(hi,1.)
        np.testing.assert_allclose(M.clopper_pearson(5,10,.05),[.18708602844739855,.8129139715526015])

    def test_exact_coverage_and_bonferroni_family_guarantee(self):
        alpha=.05/35
        intervals=[M.clopper_pearson(k,12,alpha) for k in range(13)]
        coverages=[]
        for p in (0.,.001,.02,.1,.2,.4,.5,.7,.9,.99,1.):
            covered=[lo<=p<=hi for lo,hi in intervals]
            coverage=float(np.dot(binom.pmf(np.arange(13),12,p),covered))
            self.assertGreaterEqual(coverage,1-alpha-1e-12)
            coverages.append(coverage)
        self.assertGreaterEqual(min(coverages)**35,.95-1e-12)

    def test_intervals_reject_invalid_counts(self):
        for k,n,alpha in ((-1,10,.05),(11,10,.05),(0,0,.05),(True,10,.05),(1.5,10,.05),(1,10,0)):
            with self.subTest(k=k,n=n,alpha=alpha),self.assertRaises(ValueError):
                M.clopper_pearson(k,n,alpha)

    def test_linear_compatible_regions_keep_disjoint_components(self):
        np.testing.assert_allclose(M.compatible_components([0.,1.,2.],[-2.,1.,-2.],[-1.,2.,-1.]),
                                   [[1/3,2/3],[4/3,5/3]],atol=1e-15,rtol=0)
        self.assertEqual(M.compatible_components([0.,1.],[1.,2.],[2.,3.]),[])
        self.assertEqual(M.compatible_components([0.,1.],[-1.,-1.],[1.,1.]),[[0.,1.]])

    def test_linear_regions_keep_boundary_and_point_roots(self):
        self.assertEqual(M.compatible_components([0.,1.],[0.,1.],[0.,2.]),[[0.,0.]])
        self.assertEqual(M.compatible_components([0.,1.],[-1.,0.],[-1.,0.]),[[1.,1.]])
        self.assertEqual(M.compatible_components([0.,1.,2.],[-1.,0.,1.],[-1.,0.,1.]),[[1.,1.]])
        self.assertEqual(M.intersect_components([[0.,.2],[.4,.8]],[[.1,.5],[.7,1.]]),
                         [[.1,.2],[.4,.5],[.7,.8]])

    def test_alternative_gradient_and_hessian(self):
        t=np.linspace(-1,1,7);design=np.column_stack([np.ones(7),t,t*t])
        shots=np.full(7,100.);counts=np.array([10,18,30,42,51,60,72.]);b=np.array([-.2,.9,-.1])
        loss,gradient,hessian=M.logistic_objective(b,design,counts,shots)
        eps=1e-5;fd=[];hd=[]
        for j in range(3):
            step=np.zeros(3);step[j]=eps
            plus=M.logistic_objective(b+step,design,counts,shots)
            minus=M.logistic_objective(b-step,design,counts,shots)
            fd.append((plus[0]-minus[0])/(2*eps));hd.append((plus[1]-minus[1])/(2*eps))
        np.testing.assert_allclose(gradient,fd,rtol=2e-7,atol=1e-9)
        np.testing.assert_allclose(hessian,np.array(hd).T,rtol=2e-7,atol=1e-9)
        self.assertGreater(np.linalg.eigvalsh(hessian).min(),0.)

    def test_independent_curves_reproduce_known_quadratic_logits(self):
        rates=np.linspace(.01,.07,7);coefficients=np.array([-.4,1.1,-.15]);t=np.linspace(-1,1,7)
        points=[dict(channel='x_only',steps=L,p=float(p),shots=10**8,
                     failures=int(round(10**8*expit(coefficients@[1,u,u*u]))),
                     none_failures=int(round(10**8*expit(coefficients@[1,u,u*u]))))
                for L in (3,5) for p,u in zip(rates,t)]
        fitted=M.fit_independent(points)
        self.assertTrue(fitted['optimizer_success'])
        for row in fitted['curves']:np.testing.assert_allclose(row['coefficients'],coefficients,atol=2e-7)
        self.assertLess(fitted['deviance'],1e-5)

    def test_common_model_is_exactly_nested_and_observed_deviance_decreases(self):
        rates=np.linspace(.03,.07,7)
        points=[]
        for L in (3,5,7,9,11):
            for p in rates:
                q=expit(-.4+50*(p-.05)*(L/7)**.8)
                k=int(round(3000*q))
                points.append(dict(channel='x_only',steps=L,p=float(p),shots=3000,failures=k,none_failures=k))
        common=M.fit_common(points)
        alternative=M.fit_independent(points,common)
        self.assertTrue(common['optimizer_success']);self.assertTrue(alternative['optimizer_success'])
        self.assertLessEqual(alternative['deviance'],common['deviance']+1e-7)
        np.testing.assert_allclose(M.mapped_common_probabilities(points,common),common['fitted_probabilities'],atol=1e-14)
        self.assertAlmostEqual(common['pc'],.05,delta=2e-4)

    def test_likelihood_ratio_detects_strong_size_specific_curve_distortion(self):
        points=[]
        for L,offset,slope in ((3,-1.2,.4),(5,.8,1.7),(7,-.3,.8)):
            for t in np.linspace(-1,1,7):
                k=int(round(100000*expit(offset+slope*t-.1*t*t)))
                points.append(dict(channel='x_only',steps=L,p=float(.05+.02*t),shots=100000,
                                   failures=k,none_failures=k))
        common=M.fit_common(points);alternative=M.fit_independent(points,common)
        self.assertTrue(alternative['optimizer_success'])
        self.assertLess(alternative['deviance'],.001)
        self.assertGreater(common['deviance']-alternative['deviance'],1000.)

    def test_tail_probability_reports_plus_one_and_failure_bounds(self):
        result=M.tail_summary(4,20,25)
        self.assertAlmostEqual(result['p_plus_one'],5/21)
        np.testing.assert_allclose(result['all_attempt_failure_bounds'],[5/26,10/26])
        self.assertEqual(result['numerical_failures'],5)
        self.assertLess(result['mc_binomial_interval'][0],.2)
        self.assertGreater(result['mc_binomial_interval'][1],.2)


if __name__=='__main__':unittest.main()
