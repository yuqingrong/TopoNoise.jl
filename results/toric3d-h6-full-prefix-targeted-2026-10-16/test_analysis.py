#!/usr/bin/env python3
"""Focused integrity and numerical-glue tests; no circuit simulation."""
import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import binom

HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE/'selection.json').read_text())


def load_analysis():
    path = HERE/'analyze.py'
    if not path.exists():
        return None
    spec = importlib.util.spec_from_file_location('targeted_analysis_test', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def fixture_points():
    rows = []
    for job in PLAN['jobs']:
        for size in job['sizes']:
            for p in job['rates']:
                n = job['shots']; k = n//4
                rows.append(dict(channel=job['channel'],steps=size,Ly=size,Lz=size,p=p,
                    shots=n,target_shots=n,seed=job['seed'],point_seed=f'{len(rows)+1:016x}',
                    complete=True,failures=k,none_failures=k,full_prefix_failures=k,
                    final_only_failures=2*k,full_prefix_logical_1_failures=k,
                    full_prefix_logical_2_failures=0,final_only_logical_1_failures=2*k,
                    final_only_logical_2_failures=0,joint=[n-2*k,k,0,k],source_job=job['name']))
    return rows


def completion_fixture():
    # Literal four-fit/one-profile/two-replicate fixture, independent of schedules.
    plan={'analysis':{
        'x':{'min_sizes':[5],'windows':{'full':[.0475,.0625]},'models':['quadratic'],'bootstrap_replicates':2},
        'z':{'min_sizes':[5],'windows':{'full':[.0275,.04]},'models':['quadratic'],'profile_models':['quadratic']}}}
    fit={'fit_id':'Lmin5-full-quadratic','optimizer_success':True,'pc':.034,'nu':1.2,
         'deviance':2.,'nll':100.,'bound_hits':[],'monotone_increasing':True,'deviance_p_value':.5}
    profile={'fit_id':'Lmin5-full-quadratic','minimum_agreement':True,'numerical_failures':[],
             'evaluations':[{'optimizer_success':True,'nll':100.,'lr':0.,'bound_hits':[]}],
             'regions':[{'lower':.033,'upper':.035,'lower_numerical_gap':False,'upper_numerical_gap':False,
                         'lower_bound_truncated':False,'upper_bound_truncated':False,'refinement_failures':[]}],
             'nuisance_bound_evaluation_count':0}
    channel={'fits':[fit],'profiles':[],'failures':[]}
    report={'fresh':{'x':copy.deepcopy(channel),'z':copy.deepcopy(channel)},
            'pooled_secondary':{'x':copy.deepcopy(channel),'z':copy.deepcopy(channel)}}
    report['fresh']['z']['profiles']=[profile]
    report['fresh']['x'].update(common_model=copy.deepcopy(fit),independent_size_model={'optimizer_success':True,'deviance':1.},
        bootstrap={'replicates':2,'records':[
            {'replicate':0,'common_valid':True,'alternative_valid':True,'pc':.034,'nu':1.2,'common_deviance':2.,'alternative_deviance':1.,'likelihood_ratio':1.,'bound_hits':[],'monotone':True},
            {'replicate':1,'common_valid':True,'alternative_valid':True,'pc':.034,'nu':1.2,'common_deviance':2.,'alternative_deviance':1.,'likelihood_ratio':1.,'bound_hits':[],'monotone':True}]})
    return report,plan


class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = load_analysis()

    def setUp(self):
        self.assertIsNotNone(self.a, 'Analysis module is not implemented')
        self.points = fixture_points()

    def test_failed_free_profile_point_is_retained_and_later_profiles_continue(self):
        from unittest.mock import patch
        data=self.a.select_data(self.points,PLAN,'z_only',5,'full')
        fit={'fit_id':'Lmin5-full-quadratic','model':'quadratic','min_size':5,
             'pc':.034,'nll':100.,'parameters':[.52,0.,0.,2.,0.]}
        fail_free=[True]
        def failed_at_free(data,model,pc,warm_starts=()):
            if fail_free[0] and pc==.034:
                raise ValueError('induced free-point failure')
            return {'pc':pc,'nll':100.+1e6*(pc-.034)**2,'parameters':[.52,0.,0.,2.,0.],
                    'optimizer_success':True,'optimizer_message':'synthetic optimum',
                    'failed_starts':0,'nonfinite_starts':0,'projected_gradient_infinity_norm':0.,
                    'polish_success':True,'bound_hits':[]}
        # Real reviewed profile loop and failure path; only nuisance optimization is substituted.
        with patch.object(self.a.Z,'fit_fixed_pc',side_effect=failed_at_free):
            with self.assertRaises(TypeError):
                self.a.Z.profile_pc(data,'quadratic',fit,31)
            result=self.a.profile_for_fit(data,'quadratic',fit,31)
        self.assertEqual(result['fit_id'],'Lmin5-full-quadratic')
        self.assertTrue(result['critical_helper_exception'])
        self.assertFalse(result['minimum_agreement'])
        self.assertEqual(result['numerical_failures'][0]['exception_type'],'TypeError')
        self.assertIn('NoneType',result['numerical_failures'][0]['message'])
        self.assertEqual(result['evaluations'],[])
        report,plan=completion_fixture();report['fresh']['z']['profiles']=[result]
        summary=self.a.summarize_numerics(report,plan)
        self.assertEqual(summary['status'],'complete_with_numerical_failures')
        self.assertEqual(summary['profiles']['critical_helper_exception_count'],1)
        # Retaining the first failure lets the caller continue its schedule.
        fail_free[0]=False
        with patch.object(self.a.Z,'fit_fixed_pc',side_effect=failed_at_free):
            later=self.a.profile_for_fit(dict(data,window='narrow'),'quadratic',dict(fit,fit_id='Lmin5-narrow-quadratic'),31)
        self.assertTrue(later['minimum_agreement'])
        self.assertEqual(later['window'],'narrow')

    def test_numerical_completion_requires_every_planned_result(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        report,plan=completion_fixture()
        result=self.a.summarize_numerics(report,plan)
        self.assertEqual(result['status'],'complete')
        self.assertEqual(result['fits']['expected_count'],4)
        self.assertEqual(result['profiles']['expected_count'],1)
        self.assertEqual(result['bootstrap']['expected_replicates'],2)
        report['pooled_secondary']['z']['fits']=[]
        result=self.a.summarize_numerics(report,plan)
        self.assertEqual(result['status'],'complete_with_numerical_failures')
        self.assertEqual(result['fits']['missing_count'],1)

    def test_returned_optimizer_failure_changes_completion(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        report,plan=completion_fixture()
        report['fresh']['x']['fits'][0]['optimizer_success']=False
        result=self.a.summarize_numerics(report,plan)
        self.assertEqual(result['status'],'complete_with_numerical_failures')
        self.assertEqual(result['fits']['optimizer_failure_count'],1)
        self.assertEqual(result['exception_count'],0)

    def test_profile_numerical_failures_change_completion_without_exception(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        for field,value in [('numerical_failures',[{'pc':.034,'message':'failed'}]),('minimum_agreement',False)]:
            with self.subTest(field=field):
                report,plan=completion_fixture();report['fresh']['z']['profiles'][0][field]=value
                result=self.a.summarize_numerics(report,plan)
                self.assertEqual(result['status'],'complete_with_numerical_failures')
                self.assertEqual(result['profiles']['invalid_profile_count'],1)

    def test_invalid_or_missing_bootstrap_replicate_changes_completion(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        for mode in ['common','alternative','missing','duplicate']:
            with self.subTest(mode=mode):
                report,plan=completion_fixture();records=report['fresh']['x']['bootstrap']['records']
                if mode=='missing': records.pop()
                elif mode=='duplicate': records[1]['replicate']=0
                else: records[0][mode+'_valid']=False
                result=self.a.summarize_numerics(report,plan)
                self.assertEqual(result['status'],'complete_with_numerical_failures')

    def test_statistical_flags_remain_distinct_from_numerical_failures(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        report,plan=completion_fixture();fit=report['fresh']['x']['fits'][0]
        fit.update(bound_hits=['pc'],monotone_increasing=False,deviance_p_value=.001)
        report['fresh']['z']['profiles'][0]['nuisance_bound_evaluation_count']=2
        record=report['fresh']['x']['bootstrap']['records'][0];record.update(bound_hits=['nu'],monotone=False)
        result=self.a.summarize_numerics(report,plan)
        self.assertEqual(result['status'],'complete')
        self.assertEqual(result['flags']['fit_bound_hit_count'],1)
        self.assertEqual(result['flags']['fit_nonmonotone_count'],1)
        self.assertEqual(result['flags']['fit_gof_p_lt_0p01_count'],1)
        self.assertEqual(result['flags']['profile_nuisance_bound_evaluation_count'],2)
        self.assertEqual(result['flags']['bootstrap_bound_hit_count'],1)
        self.assertEqual(result['flags']['bootstrap_nonmonotone_count'],1)

    def test_nonfinite_numeric_results_and_profile_gaps_change_completion(self):
        self.assertTrue(callable(getattr(self.a,'summarize_numerics',None)), 'numerical completion aggregator missing')
        report,plan=completion_fixture();report['fresh']['x']['fits'][0]['pc']=float('nan')
        result=self.a.summarize_numerics(report,plan)
        self.assertEqual(result['status'],'complete_with_numerical_failures')
        self.assertEqual(result['fits']['nonfinite_fit_count'],1)
        report,plan=completion_fixture();report['fresh']['z']['profiles'][0]['regions'][0]['lower_numerical_gap']=True
        self.assertEqual(self.a.summarize_numerics(report,plan)['status'],'complete_with_numerical_failures')

    def test_frozen_grid_totals_and_unequal_counts(self):
        result = self.a.validate_fresh_points(self.points, PLAN)
        self.assertEqual(result['records'], 204000)
        self.assertEqual(result['by_channel'], {'x_only':106000,'z_only':98000})
        self.assertEqual(result['cells'], 77)
        data = self.a.select_data(self.points, PLAN, 'x_only', 5, 'full')
        self.assertEqual(data['shots'].sum(), 106000)
        self.assertEqual(set(data['shots']), {2000,4000,5000})

    def test_missing_cell_rejected(self):
        with self.assertRaisesRegex(ValueError, 'grid'):
            self.a.validate_fresh_points(self.points[:-1], PLAN)

    def test_incomplete_count_rejected(self):
        self.points[0]['shots'] -= 1
        with self.assertRaisesRegex(ValueError, 'shot|count|complete'):
            self.a.validate_fresh_points(self.points, PLAN)

    def test_incomplete_flag_rejected(self):
        self.points[0]['complete'] = False
        with self.assertRaisesRegex(ValueError, 'complete'):
            self.a.validate_fresh_points(self.points, PLAN)

    def test_duplicate_cell_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate|grid'):
            self.a.validate_fresh_points(self.points+[self.points[0]], PLAN)

    def test_duplicate_stream_rejected_across_different_cells(self):
        self.points[1]['point_seed'] = self.points[0]['point_seed']
        with self.assertRaisesRegex(ValueError, 'stream'):
            self.a.validate_fresh_points(self.points, PLAN)

    def test_old_stream_collision_rejected_even_outside_pool(self):
        old = dict(self.points[0],steps=3,Ly=3,Lz=3,p=.001,source_job='old')
        with self.assertRaisesRegex(ValueError, 'stream'):
            self.a.assert_independent(self.points,[old])

    def test_excluded_timing_stream_collision_rejected(self):
        with self.assertRaisesRegex(ValueError, 'stream'):
            self.a.assert_independent(self.points,[],[self.points[0]['point_seed']])

    def test_wrong_job_seed_rejected(self):
        self.points[0]['seed'] += 1
        with self.assertRaisesRegex(ValueError, 'seed'):
            self.a.validate_fresh_points(self.points, PLAN)

    def test_joint_histogram_inconsistency_rejected(self):
        self.points[0]['joint'][0] -= 1
        with self.assertRaisesRegex(ValueError, 'joint|count'):
            self.a.validate_fresh_points(self.points, PLAN)

    def test_pooled_counts_add_only_matching_fresh_cells(self):
        old = [dict(self.points[0],shots=1000,target_shots=1000,failures=300,
                    full_prefix_failures=300,none_failures=300,point_seed='fffffffffffffffe',source_job='old'),
               dict(self.points[0],steps=3,Ly=3,Lz=3,shots=1000,point_seed='ffffffffffffffff',source_job='old')]
        pooled, info = self.a.pool_points(self.points, old)
        row = next(r for r in pooled if self.a.cell_key(r)==self.a.cell_key(self.points[0]))
        self.assertEqual(row['shots'],3000)
        self.assertEqual(row['failures'],800)
        self.assertEqual(info['old_records_added'],1000)
        self.assertEqual(sum(r['shots'] for r in pooled),205000)
        self.assertEqual(len(pooled),77)

    def test_window_mapping_is_channel_local(self):
        x = self.a.select_data(self.points,PLAN,'x_only',7,'narrow')
        z = self.a.select_data(self.points,PLAN,'z_only',9,'narrow')
        x2 = self.a.select_data(self.points,PLAN,'x_only',7,'narrow')
        self.assertEqual((x['p_min'],x['p_max'],len(x['p'])),(.05,.06,15))
        self.assertEqual((z['p_min'],z['p_max'],len(z['p'])),(.03,.0375,25))
        np.testing.assert_array_equal(x['p'],x2['p'])
        self.assertEqual(self.a.Z.WINDOWS, {'full':(.025,.04),'narrow':(.0275,.0375)})

    def test_model_and_profile_schedules(self):
        self.assertEqual(len(self.a.fit_schedule(PLAN,'x_only')),24)
        self.assertEqual(len(self.a.fit_schedule(PLAN,'z_only')),36)
        self.assertEqual(len(self.a.profile_schedule(PLAN)),12)
        self.assertIn((9,'narrow','corr_w1'),self.a.profile_schedule(PLAN))

    def test_anchor_uses_four_cells_and_correct_direction(self):
        for row in self.points:
            if row['channel']=='x_only' and row['steps'] in [5,11] and row['p'] in [.05,.06]:
                row['failures']=int(row['shots']*({(5,.05):.25,(11,.05):.10,(5,.06):.30,(11,.06):.50}[row['steps'],row['p']]))
        result = self.a.anchor_assessment(self.points,PLAN)
        self.assertEqual(result['number_cells'],4)
        self.assertAlmostEqual(result['cell_alpha'],.0125)
        self.assertTrue(result['success'])
        self.assertLess(result['contrasts'][0]['ci_high'],0)
        self.assertGreater(result['contrasts'][1]['ci_low'],0)
        self.assertEqual([c['shots'] for c in result['cells']],[5000]*4)

    def test_anchor_equal_probabilities_do_not_resolve_crossing(self):
        self.assertFalse(self.a.anchor_assessment(self.points,PLAN)['success'])

    def test_exact_anchor_cell_coverage_small_n(self):
        # Exact enumeration, not Monte Carlo: CP cell coverage >= .9875.
        for p in [.001,.1,.5,.9,.999]:
            coverage=0.
            for k in range(21):
                lo,hi=self.a.X.clopper_pearson(k,20,.0125)
                coverage += binom.pmf(k,20,p)*(lo<=p<=hi)
            self.assertGreaterEqual(coverage, .9875-1e-12)

    def test_secondary_bands_describe_actual_grid(self):
        result=self.a.secondary_crossing([r for r in self.points if r['channel']=='x_only'])
        self.assertEqual(result['number_cells'],28)
        self.assertIn('28',result['method'])
        self.assertEqual(len(result['pairs']),6)
        self.assertAlmostEqual(result['cell_alpha'],.05/28)

    def test_old_shard_integrity_and_source_bound_loading(self):
        old=HERE.parent/'toric3d-h6-full-prefix-2026-10-15'
        rows, provenance=self.a.load_old_production(PLAN,old)
        self.assertEqual(len(rows),70)
        self.assertEqual(sum(r['shots'] for r in rows),70000)
        self.assertEqual(len(provenance['shards']),6)
        self.assertFalse(any('pilot' in r['source_job'] for r in rows))

    def test_selection_digest_rejects_changed_plan(self):
        self.assertTrue(callable(getattr(self.a,'read_plan',None)), 'frozen plan digest reader missing')
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            (directory/'selection.json').write_bytes((HERE/'selection.json').read_bytes())
            (directory/'selection.sha256').write_bytes((HERE/'selection.sha256').read_bytes())
            self.assertEqual(self.a.read_plan(directory)['total_new_fault_records'],204000)
            with (directory/'selection.json').open('a') as stream: stream.write(' ')
            with self.assertRaisesRegex(ValueError,'selection.*hash|digest'):
                self.a.read_plan(directory)

    def test_excluded_benchmark_loads_with_actual_runner_schema(self):
        import tomllib
        directory=HERE/'benchmark-z'
        metadata=tomllib.loads((directory/'metadata.toml').read_text())
        reference=tomllib.loads((HERE.parent/'toric3d-h6-full-prefix-2026-10-15'/'production-z'/'metadata.toml').read_text())
        job=dict(metadata['config'],name='benchmark-z')
        rows,provenance=self.a.read_shard(directory,job,PLAN['simulation_source_sha256'],reference['protocol'])
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['shots'],20)
        self.assertEqual(rows[0]['steps'],17)
        self.assertEqual(provenance['base_seed'],202610166002)

    def test_raw_checksum_tampering_rejected(self):
        import shutil,tomllib
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)/'shard';shutil.copytree(HERE/'benchmark-z',directory)
            metadata=tomllib.loads((directory/'metadata.toml').read_text())
            with (directory/'raw.csv').open('a') as stream: stream.write('\n')
            with self.assertRaisesRegex(ValueError,'raw checksum'):
                self.a.read_shard(directory,dict(metadata['config'],name='benchmark-z'))

    def test_narrow_profile_preserves_actual_window_and_contains_truth(self):
        self.assertTrue(callable(getattr(self.a,'profile_for_fit',None)), 'window-safe profile wrapper missing')
        points=copy.deepcopy(self.points)
        for row in points:
            if row['channel']=='z_only':
                u=(row['p']-.034)/.0125*(row['steps']/10)**(1/1.2)
                row['shots']=100000;row['failures']=round(100000/(1+np.exp(-(.1+2*u-.3*u*u))))
        data=self.a.select_data(points,PLAN,'z_only',5,'narrow')
        fit=self.a.Z.fit_model(data,'quadratic');fit['fit_id']='Lmin5-narrow-quadratic'
        profile=self.a.profile_for_fit(data,'quadratic',fit,31)
        self.assertEqual(profile['window'],'narrow')
        self.assertEqual(profile['pc_bounds'],[.03,.0375])
        self.assertTrue(profile['minimum_agreement'])
        self.assertTrue(any(r['lower']<.034<r['upper'] for r in profile['regions']))

    def test_bootstrap_uses_frozen_seed_unequal_counts_and_isolated_settings(self):
        self.assertTrue(callable(getattr(self.a,'run_bootstrap',None)), 'isolated bootstrap wrapper missing')
        points=[copy.deepcopy(r) for r in self.points if r['channel']=='x_only']
        for row in points:
            u=(row['p']-.054)/.015*(row['steps']/8)**(1/1.4)
            row['failures']=round(row['shots']/(1+np.exp(-(-1.1+1.2*u-.2*u*u))))
            row['none_failures']=row['failures']
        plan=copy.deepcopy(PLAN);plan['analysis']['x']['bootstrap_replicates']=2
        original=(self.a.X.B,self.a.X.SEED)
        result=self.a.run_bootstrap(points,plan)
        boot=result['bootstrap'];common=result['common_model']
        expected=np.random.default_rng(202610166901).binomial(np.array([r['shots'] for r in points]),common['fitted_probabilities'],size=(2,len(points)))
        self.assertEqual(boot['replicates'],2)
        np.testing.assert_array_equal([r['simulated_failures'] for r in boot['records']],expected)
        self.assertEqual((self.a.X.B,self.a.X.SEED),original)
        self.assertEqual(boot['goodness_of_fit']['valid_replicates'],2)

    def test_numerical_sources_match_previous_verified_reports(self):
        source_hashes=self.a.verify_numerical_sources(PLAN)
        self.assertIn('results/toric3d-h6-full-prefix-targeted-2026-10-16/analyze.py',source_hashes)
        self.assertIn('examples/fit_toric3d_h6_threshold.py',source_hashes)

    def test_excluded_timing_and_pilot_streams_are_auditable(self):
        seeds,records=self.a.collect_excluded_streams(HERE,self.a.OLD,PLAN)
        self.assertEqual(len(records),5)
        self.assertGreater(len(set(seeds)),2)
        self.assertEqual(sum('benchmark-' in r['path'] for r in records),2)
        self.assertEqual(sum('/pilot' in r['path'] for r in records),3)
        self.assertTrue(all(len(r['sha256'])==64 for r in records))

    def test_synthetic_unequal_shot_fit_recovers_probability_curve(self):
        points=copy.deepcopy(self.points)
        for row in points:
            if row['channel']=='z_only':
                u=(row['p']-.034)/.0125*(row['steps']/10)**(1/1.2)
                q=1/(1+np.exp(-(.1+2*u-.3*u*u)))
                row['shots']=100000+10000*row['steps'];row['failures']=round(row['shots']*q);row['full_prefix_failures']=row['failures']
        data=self.a.select_data(points,PLAN,'z_only',5,'full')
        fit=self.a.Z.fit_model(data,'quadratic')
        self.assertTrue(fit['optimizer_success'])
        self.assertAlmostEqual(fit['pc'],.034,delta=2e-6)
        self.assertAlmostEqual(fit['nu'],1.2,delta=.002)
        self.assertLess(max(abs(np.array(fit['fitted_probabilities'])-data['failures']/data['shots'])),2e-5)


if __name__=='__main__': unittest.main(verbosity=2)
