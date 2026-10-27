import copy
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import tempfile
import time
import unittest

PATH = Path(__file__).with_name('analyze.py')
M = None
if PATH.exists():
    spec = importlib.util.spec_from_file_location('axis_analysis', PATH)
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)


def selection():
    channels = {ch: dict(sizes=[3], rates=[.005, .01]) for ch in ('x_only', 'z_only')}
    return dict(stage='direction_comparison', frozen_before_sampling=True,
        frozen_at_unix=100., shots=1000, batch_size=25,
        arms=['final_only', 'one_step', 'two_slice'], channels=channels,
        common_z_seed_repair=True, pilot_excluded=True,
        total_fault_records=4000, total_protocol_outcomes=12000,
        shards={f'{ch}-L3': dict(channel=ch, size=3, rates=[.005,.01], seed=100+i,
                  point_seeds_hex=[f'{10+2*i+j:016x}' for j in range(2)])
                for i,ch in enumerate(channels)},
        protocols={ch: dict(channel=ch, arm_order=['final_only','one_step','two_slice']) for ch in channels})


def row(sel, channel='x_only', p=.005, shots=1000):
    shard = sel['shards'][f'{channel}-L3']
    result = dict(channel=channel, steps=3, Ly=3, Lz=3, p=p, shots=shots,
        target_shots=1000, seed=shard['seed'],
        point_seed=shard['point_seeds_hex'][shard['rates'].index(p)], complete=shots==1000,
        elapsed_seconds=1., source_shard=f'{channel}-L3')
    for arm in ('final_only','one_step','two_slice'):
        result.update({arm+'_logical_1_failures':0,arm+'_logical_2_failures':0,arm+'_failures':0})
    result.update({f'joint_{mask:02}':shots if mask==0 else 0 for mask in range(8)})
    return result


def toml_value(value):
    if isinstance(value,dict):
        return '{'+', '.join(json.dumps(k)+' = '+toml_value(v) for k,v in sorted(value.items()))+'}'
    if isinstance(value,list):
        return '['+', '.join(toml_value(v) for v in value)+']'
    return json.dumps(value)


def toml_text(data):
    return ''.join(k+' = '+toml_value(v)+'\n' for k,v in sorted(data.items()))


def write_fixture(folder, sel, name, *, partial=False):
    shard=sel['shards'][name]; channel=shard['channel']
    rows=[row(sel,channel,p,25 if partial and i==1 else 1000) for i,p in enumerate(shard['rates'])]
    runtime=dict(julia_threads=1,blas_threads=1,python_executable='/fixture/python',
        julia_version='test',packages={k:'test' for k in ('numpy','scipy','ldpc','pymatching')},
        environment={**{k:'1' for k in ('JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
            'JULIA_CONDAPKG_BACKEND':'Null','JULIA_PYTHONCALL_EXE':'/fixture/python'})
    config=dict(channel=channel,sizes=[3],rates=shard['rates'],shots=1000,seed=shard['seed'],
        batch_size=25,dimensions=[[3,3,3]])
    states=[dict(completed_shots=r['shots'],elapsed_seconds=1.,logical_1_counts=[0,0,0],
        logical_2_counts=[0,0,0],failures=[0,0,0],joint_counts=[r['shots'],0,0,0,0,0,0,0]) for r in rows]
    cp=dict(format_version=1,config=config,runtime=runtime,source_hashes=sel.get('source_hashes',{}),
        protocol=sel['protocols'][channel],point_seeds_hex=shard['point_seeds_hex'],states=states,
        generation=2,started_at_unix=101.,updated_at_unix=102.,
        decoder_build_seconds_by_size={'3':1.},decoder_build_count_by_size={'3':2})
    payload=toml_text(cp)
    checkpoint='payload_sha256 = "'+hashlib.sha256(payload.encode()).hexdigest()+'"\n'+payload
    stream=io.StringIO(); writer=csv.DictWriter(stream,fieldnames=M.COLUMNS); writer.writeheader(); writer.writerows(rows)
    raw=stream.getvalue()
    meta={k:v for k,v in cp.items() if k!='states'}
    meta.update(expected_points=2,completed_points=sum(r['complete'] for r in rows),
        completed_shots=sum(r['shots'] for r in rows),complete=all(r['complete'] for r in rows),
        checkpoint_sha256=hashlib.sha256(checkpoint.encode()).hexdigest(),raw_sha256=hashlib.sha256(raw.encode()).hexdigest())
    folder.mkdir(parents=True)
    (folder/'checkpoint.toml').write_text(checkpoint); (folder/'metadata.toml').write_text(toml_text(meta))
    (folder/'raw.csv').write_text(raw)
    return cp,rows


class ImplementationTest(unittest.TestCase):
    def test_analysis_implementation_exists(self):
        self.assertTrue(PATH.is_file(), 'direction analysis is not implemented yet')


@unittest.skipUnless(M is not None, 'implementation not written yet')
class AnalysisTests(unittest.TestCase):
    def test_invalid_counts_are_rejected(self):
        sel=selection(); good=row(sel)
        self.assertEqual(M.validate_row(good,sel)['shots'],1000)
        for field,value in [('joint_00',999),('final_only_failures',1),('shots',True),
                            ('one_step_logical_2_failures',-1),('elapsed_seconds',math.nan)]:
            bad=dict(good); bad[field]=value
            with self.assertRaises(ValueError,msg=field): M.validate_row(bad,sel)

    def test_wrong_master_and_point_seeds_are_rejected(self):
        sel=selection()
        for field,value in [('seed',999),('point_seed','ffffffffffffffff')]:
            bad=row(sel); bad[field]=value
            with self.assertRaisesRegex(ValueError,'seed'): M.validate_row(bad,sel)
        bad=copy.deepcopy(sel)
        bad['shards']['z_only-L3']['point_seeds_hex'][0]=bad['shards']['x_only-L3']['point_seeds_hex'][0]
        with self.assertRaisesRegex(ValueError,'seed'): M.validate_selection(bad)
        bad=copy.deepcopy(sel); bad['shards']['z_only-L3']['rates']=[.005]
        with self.assertRaises(ValueError): M.validate_selection(bad)

    def test_independent_direction_difference_uses_sum_of_binomial_variances(self):
        answer=M.independent_difference(200,1000,100,1000)
        self.assertAlmostEqual(answer['difference'],.1)
        self.assertAlmostEqual(answer['standard_error'],math.sqrt(.2*.8/1000+.1*.9/1000))
        self.assertAlmostEqual(answer['ci_low'],0.06894934015909647)
        self.assertAlmostEqual(answer['ci_high'],0.13104641439407827)
        self.assertEqual(answer['sampling'],'independent across directions')
        self.assertIn('Newcombe',answer['interval_method'])
        with self.assertRaises(ValueError): M.independent_difference(1,0,1,10)

    def test_independent_score_interval_handles_boundaries_unequal_n_and_symmetry(self):
        zero=M.independent_difference(0,1000,0,1000)
        self.assertEqual(zero['standard_error'],0.)
        self.assertAlmostEqual(zero['ci_low'],-0.003826758485555124)
        self.assertAlmostEqual(zero['ci_high'],0.003826758485555124)
        self.assertGreater(zero['ci_high']-zero['ci_low'],0.)
        full=M.independent_difference(1000,1000,1000,1000)
        self.assertAlmostEqual(full['ci_low'],zero['ci_low'])
        self.assertAlmostEqual(full['ci_high'],zero['ci_high'])
        unequal=M.independent_difference(100,200,50,100)
        self.assertAlmostEqual(unequal['standard_error'],math.sqrt(.25/200+.25/100))
        self.assertAlmostEqual(unequal['ci_low'],-0.11815120035811344)
        self.assertAlmostEqual(unequal['ci_high'],0.1181512003581135)
        forward=M.independent_difference(20,200,60,1000)
        reverse=M.independent_difference(60,1000,20,200)
        self.assertAlmostEqual(forward['difference'],-reverse['difference'])
        self.assertAlmostEqual(forward['ci_low'],-reverse['ci_high'])
        self.assertAlmostEqual(forward['ci_high'],-reverse['ci_low'])

    def test_within_axis_contrasts_use_discordant_joint_bins(self):
        r=row(selection()); r['shots']=100
        r.update({f'joint_{i:02}':v for i,v in enumerate([50,10,20,5,10,1,2,2])})
        contrast=M.paired_contrasts(r)[0]
        self.assertEqual(contrast['reference'],'final_only')
        self.assertEqual(contrast['comparison'],'one_step')
        self.assertEqual(contrast['comparison_only'],22)
        self.assertEqual(contrast['reference_only'],11)
        self.assertAlmostEqual(contrast['difference'],.11)

    def test_grouped_shard_rows_and_frozen_seed_order_match_checkpoint(self):
        sel=selection(); M.validate_selection(sel)
        with tempfile.TemporaryDirectory() as d:
            study=Path(d)
            for name in sel['shards']: write_fixture(study/'archive/runs'/name,sel,name)
            rows,provenance,runtime,skipped=M.load_inputs(study,sel,{},progress=False)
            self.assertEqual(len(rows),4); self.assertEqual(len(provenance),2); self.assertEqual(skipped,[])
            folder=study/'archive/runs/x_only-L3'
            cp=(folder/'checkpoint.toml').read_text()
            damaged=cp.replace('"completed_shots" = 1000','"completed_shots" = 999',1)
            self.assertNotEqual(damaged,cp)
            (folder/'checkpoint.toml').write_text(damaged)
            with self.assertRaisesRegex(ValueError,'checksum'): M.load_inputs(study,sel,{},progress=False)

    def test_incomplete_final_is_rejected_but_progress_is_explicit(self):
        sel=selection(); rows=[row(sel)]
        with self.assertRaisesRegex(ValueError,'incomplete'): M.require_complete(rows,sel,progress=False)
        self.assertFalse(M.require_complete(rows,sel,progress=True))
        full=[row(sel,ch,p) for ch in sel['channels'] for p in sel['channels'][ch]['rates']]
        self.assertTrue(M.require_complete(full,sel,progress=False))

    def test_recent_unsynchronized_mirrors_are_skipped_only_for_progress(self):
        sel=selection()
        with tempfile.TemporaryDirectory() as d:
            study=Path(d); folder=study/'archive/runs/x_only-L3'
            write_fixture(folder,sel,'x_only-L3')
            (folder/'raw.csv').write_text('unsynchronized mirror')
            rows,_,_,skipped=M.load_inputs(study,sel,{},progress=True)
            self.assertEqual(rows,[]); self.assertEqual(skipped,['x_only-L3'])
            with self.assertRaises((M.ActiveSnapshot,ValueError)): M.load_inputs(study,sel,{},progress=False)

    def test_reference_hash_and_independent_rng_streams_are_enforced(self):
        sel=selection()
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d); ref_rows=[]
            for ch in sel['channels']:
                for p in sel['channels'][ch]['rates']:
                    r=row(sel,ch,p); r['seed']+=500; r['point_seed']=f'{500+len(ref_rows):016x}'; ref_rows.append(r)
            stream=io.StringIO(); writer=csv.DictWriter(stream,fieldnames=M.COLUMNS); writer.writeheader(); writer.writerows(ref_rows)
            csv_bytes=stream.getvalue().encode(); (folder/'reference.csv').write_bytes(csv_bytes)
            report={'complete':True,'output_sha256':{'data/matched-results.csv':M.digest(csv_bytes)}}
            (folder/'reference.json').write_text(json.dumps(report))
            sel['reference']=dict(csv=str(folder/'reference.csv'),csv_sha256=M.digest(csv_bytes),
                report=str(folder/'reference.json'),report_sha256=M.digest((folder/'reference.json').read_bytes()))
            self.assertEqual(len(M.load_reference(sel)),4)
            sel['shards']['x_only-L3']['point_seeds_hex'][0]=ref_rows[0]['point_seed']
            with self.assertRaisesRegex(ValueError,'independent'): M.load_reference(sel)
            sel['reference']['csv_sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'hash'): M.load_reference(sel)

    def test_frozen_sources_reject_changes_and_escaping_auxiliary_paths(self):
        with tempfile.TemporaryDirectory() as d:
            study=Path(d); archive=study/'archive'; (archive/'runtime/src').mkdir(parents=True)
            names=['runtime/Project.toml','runtime/Manifest.toml','runtime/src/main.jl',
                   'run_study.jl','AxisH6.jl','axis_geometry.jl','axis_feedback.jl']
            auxiliary=['analyze.py','supervise.py','paper-no-grid.mplstyle','test_extra.py']
            for name in names+auxiliary: (archive/name).write_text('frozen fixture\n')
            sel=dict(source_hashes={n:M.digest((archive/n).read_bytes()) for n in names},
                     auxiliary_source_hashes={n:M.digest((archive/n).read_bytes()) for n in auxiliary})
            self.assertEqual(M.verify_frozen_sources(study,sel),sel['source_hashes'])
            (archive/'runtime/src/main.jl').write_text('changed\n')
            with self.assertRaisesRegex(ValueError,'hash'): M.verify_frozen_sources(study,sel)
            (archive/'runtime/src/main.jl').write_text('frozen fixture\n')
            (study/'outside.py').write_text('outside\n')
            sel['auxiliary_source_hashes']['../outside.py']=M.digest((study/'outside.py').read_bytes())
            with self.assertRaisesRegex(ValueError,'hash'): M.verify_frozen_sources(study,sel)

    def test_report_keeps_independent_and_paired_comparisons_separate_without_fits(self):
        sel=selection(); sel.update(source_hashes={},auxiliary_source_hashes={},reference={})
        reference=[row(sel,ch,p) for ch in sel['channels'] for p in sel['channels'][ch]['rates']]
        report=M.make_report([row(sel,shots=25)],reference,sel,{},None,'fixture',skipped=['z_only-L3'])
        self.assertFalse(report['complete'])
        self.assertEqual(report['completed_fault_records'],25)
        self.assertEqual(report['completed_protocol_outcomes'],75)
        self.assertEqual(report['fits'],{})
        self.assertEqual(len(report['paired_within_axis_contrasts']),3)
        self.assertEqual(len(report['independent_axis_minus_diagonal']),3)
        self.assertEqual(report['skipped_updating_shards'],['z_only-L3'])


if __name__=='__main__': unittest.main()
