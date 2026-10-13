import copy, importlib.util, json, math, tempfile, unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('h5analysis',ROOT/'examples/analyze_toric3d_vertex_projector.py')
h5=importlib.util.module_from_spec(spec);spec.loader.exec_module(h5)

def point(L,p,rates,channel='x_only',shots=100000):
    counts=[0]*16
    ks=[round(shots*r) for r in rates]
    for left,right in zip(sorted(set([0,shots]+ks))[:-1],sorted(set([0,shots]+ks))[1:]):
        mask=sum(1<<j for j,k in enumerate(ks) if left<k);counts[mask]+=right-left
    row=dict(T=L,Lu=L,Lv=L,steps=L,p=p,channel=channel,metric='logical_x_union' if channel=='x_only' else 'logical_z_union',shots=shots,decoder_build_seconds=.1*L,sampling_seconds=0.,replay_seconds=shots*(.003 if channel=='x_only' else .005)*L**3,invalid_recoveries=0,complete=True)
    row.update({f'joint_mask_{i}':n for i,n in enumerate(counts)})
    row.update({f'{pr}_failures':k for pr,k in zip(h5.PROTOCOLS,ks)})
    return h5.normalize_row(row)

def scout():
    return [point(L,p,[1/(1+math.exp(-100*(p-.03)*L**.8))]*4,ch) for ch in h5.CHANNELS for L in (3,5,7) for p in (.02,.025,.03,.035,.04)]
META=dict(seed=1234,shots=500,decoder_settings={'max_iter':50},source_hashes={'src':'a'*64},protocols=['none','one_step','window_complete','full_prefix'],stage='scout',format_version=1,complete=True,invalid_recoveries=0,extension_round=0)
INPUTS=[dict(csv_sha256='a'*64,metadata_sha256='b'*64,path_ids=['a.csv','a.toml'],stage='scout',extension_round=0,parent_selection_sha256='')]
class H5FitTests(unittest.TestCase):
    def test_reads_four_arm_joint_counts_and_six_contrasts(self):
        r=point(3,.01,[.1,.2,.3,.4]);r.update({f'joint_mask_{i}':i+1 for i in range(16)});r['shots']=136
        for pr,k in zip(h5.PROTOCOLS,(72,76,84,100)):r[f'{pr}_failures']=k
        r=h5.normalize_row(r)
        self.assertEqual(list(r['marginals'].values()),[72,76,84,100])
        self.assertEqual(len(r['pairs']),6)
        self.assertEqual(r['pairs']['none_vs_one_step']['cells'],[28,32,36,40])
        r['none_failures']=1
        with self.assertRaises(ValueError):h5.normalize_row(r)
    def test_scout_selection_is_frozen_before_production(self):
        a=h5.select_grids(scout(),META,INPUTS)
        self.assertEqual(a['status'],'ready')
        self.assertEqual(a['channels']['x_only']['candidate_centers']['none'],[.03])
        self.assertEqual(a['channels']['x_only']['rates'],[.024,.027,.03,.033,.036])
        with self.assertRaises(ValueError):h5.select_grids(scout(),dict(META,stage='production'),INPUTS)
    def test_scout_shards_merge_exact_points(self):
        rows=scout();self.assertEqual(len(h5.merge_points([rows,copy.deepcopy(rows)])),len(rows))
        bad=copy.deepcopy(rows);bad[0]['shots']+=1
        with self.assertRaises(ValueError):h5.merge_points([rows,bad])
    def test_multiple_crossings_fund_every_candidate_window(self):
        rows=[point(L,p,[.5+(L-5)*v]*4,ch) for ch in h5.CHANNELS for L in (3,5,7) for p,v in zip((.01,.02,.04,.06,.08),(-.05,.05,.05,-.05,.05))]
        a=h5.select_grids(rows,META,INPUTS)
        for ch in h5.CHANNELS:
            v=a['channels'][ch];self.assertGreaterEqual(len(v['candidate_centers']['none']),2)
            for pc,w in zip(v['candidate_centers']['none'],v['fit_windows']['none']):
                expected=[float(format(pc*k,'.8g')) for k in (.8,.9,1,1.1,1.2)]
                self.assertTrue(set(expected)<=set(v['rates']));self.assertEqual(w,[min(expected),max(expected)])
    def test_empty_saturation_extends_then_terminates_at_bounds(self):
        rows=[point(L,p,[0.]*4,ch) for ch in h5.CHANNELS for L in (3,5,7) for p in (.01,.02)]
        a=h5.select_grids(rows,META,INPUTS);self.assertEqual(a['status'],'needs_extended_scout')
        self.assertEqual(a['extension_grids']['x_only']['rates'],[.03,.04])
        a=h5.select_grids(rows,META,INPUTS,extension_round=3)
        self.assertEqual(a['status'],'ready');self.assertTrue(a['channels']['x_only']['flags']['none']['threshold_ineligible'])
    def test_mixed_arm_extension_status_uses_active_union(self):
        rows=[point(L,p,[0.,1.,.4,.4],ch) for ch in h5.CHANNELS for L in (3,5,7) for p in (.1,.2,.3,.4,.49)]
        a=h5.select_grids(rows,META,INPUTS)
        self.assertEqual(a['status'],'needs_extended_scout')
        self.assertEqual(a['channels']['x_only']['extension_states']['none'],'exhausted_diagnostic')
        self.assertEqual(a['channels']['x_only']['extension_states']['one_step'],'requesting')
        self.assertEqual(a['extension_grids']['x_only']['rates'],[.01,.033333333])
    def test_extension_rounds_require_cumulative_parent_chain(self):
        prev=h5.select_grids(scout(),META,INPUTS)
        extension=dict(INPUTS[0],csv_sha256='c'*64,stage='extend-scout',extension_round=1,parent_selection_sha256='d'*64)
        h5.validate_ancestry(INPUTS+[extension],[META,dict(META,stage='extend-scout',extension_round=1,parent_selection_sha256='d'*64,parent_scout_inputs=INPUTS)],prev,'d'*64)
        with self.assertRaises(ValueError):h5.validate_ancestry([extension],[META],prev,'d'*64)
        with self.assertRaises(ValueError):h5.validate_ancestry(INPUTS+[dict(extension,parent_selection_sha256='e'*64)],[META,META],prev,'d'*64)
    def test_timing_rule_selects_expected_channel_sizes(self):
        a=h5.select_grids(scout(),META,INPUTS)
        self.assertEqual([d[0] for d in a['channels']['x_only']['dimensions_T_Lu_Lv']],[3,5,7,9,11])
        self.assertEqual([d[0] for d in a['channels']['z_only']['dimensions_T_Lu_Lv']],[3,5,7,9])
    def test_all_threshold_guards_withhold_threshold_claims(self):
        rows=[point(L,p,[.2]*4) for L in (3,5,7,9) for p in (.024,.027,.03,.033,.036)]
        fit=h5.fit_window(rows,'x_only','none',[.024,.036],bootstrap=0)
        self.assertIsNone(fit['pc']);self.assertTrue(fit['diagnostic'])
        for change in ({'optimizer_success':False},{'monotone_increasing':False},{'bound_hits':['pc']},{'deviance_p_value':0.}):
            c=dict(optimizer_success=True,monotone_increasing=True,bound_hits=[],pc=.03,deviance_p_value=.5);c.update(change)
            self.assertTrue(h5.engine._candidate_guard_issues(c,p_min=.02,p_max=.04,goodness_alpha=.01))
    def test_seeded_bootstrap_is_reproducible_per_channel_and_protocol(self):
        rows=[point(L,p,[1/(1+math.exp(-100*(p-.03)*L**.8))]*4) for L in (3,5,7,9) for p in (.024,.027,.03,.033,.036)]
        a=h5.fit_window(rows,'x_only','none',[.024,.036],bootstrap=4,seed=1234)
        b=h5.fit_window(rows,'x_only','none',[.024,.036],bootstrap=4,seed=1234)
        self.assertEqual(a['bootstrap'],b['bootstrap'])
        self.assertEqual(a['status'],'ok');self.assertIsNotNone(a['bootstrap']['pc_interval'])
    def test_reports_hashes_sensitivities_and_diagnostic_candidates(self):
        rows=[point(L,p,[.2]*4,ch) for ch in h5.CHANNELS for L in (3,5,7,9) for p in (.024,.027,.03,.033,.036)]
        selection=h5.select_grids(scout(),META,INPUTS)
        report=h5.analyze(rows,selection,bootstrap=0,seed=1234)
        self.assertEqual(len(report['fits']),8)
        self.assertTrue(all(f['pc'] is None for f in report['fits']))
        self.assertTrue(report['sensitivity'])
        with tempfile.TemporaryDirectory() as d:
            h5.write_report(report,rows,Path(d),hashes={'selection':'a'*64})
            self.assertTrue((Path(d)/'report.md').is_file())
            self.assertIn('diagnostic', (Path(d)/'report.md').read_text().lower())
            self.assertEqual(json.loads((Path(d)/'fits.json').read_text())['hashes']['selection'],'a'*64)
    def test_high_probability_windows_preserve_five_rate_rule(self):
        rows=[point(L,p,[.4]*4,ch) for ch in h5.CHANNELS for L in (3,5,7) for p in (.42,.44,.46,.48,.49)]
        selected=h5.select_grids(rows,META,INPUTS)
        self.assertEqual(selected['channels']['x_only']['candidate_centers']['none'],[.44])
        self.assertEqual(selected['channels']['x_only']['rates'],[.352,.396,.44,.484,.528])
        f=h5.fit_window(rows,'x_only','none',[.352,.528],bootstrap=0)
        self.assertNotIn('invalid fit window',f['diagnostic'])
        crossing=[point(L,p,[.5+(L-5)*(p-.43)]*4,ch) for ch in h5.CHANNELS for L in (3,5,7) for p in (.36,.39,.42,.45,.48)]
        selected=h5.select_grids(crossing,META,INPUTS)
        self.assertEqual(selected['channels']['x_only']['candidate_centers']['none'],[.43])
        self.assertEqual(selected['channels']['x_only']['fit_windows']['none'],[[.344,.516]])
    def test_report_refuses_a_csv_with_the_wrong_saved_hash(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'data.csv').write_text('changed data\n')
            (d/'fits.json').write_text(json.dumps({'fits':[],'sensitivity':[],'hashes':{'raw_csv':'a'*64}}))
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                h5.main(['report','--fits',str(d/'fits.json'),'--csv',str(d/'data.csv'),'--output',str(d/'out')])
            self.assertFalse((d/'out').exists())
if __name__=='__main__':unittest.main()
