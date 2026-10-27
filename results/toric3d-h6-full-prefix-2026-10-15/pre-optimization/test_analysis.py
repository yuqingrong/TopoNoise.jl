import importlib.util
import importlib.metadata
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest

PATH = Path(__file__).with_name('analyze.py')
M = None
if PATH.exists():
    spec = importlib.util.spec_from_file_location('full_prefix_analysis', PATH)
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)


def row(channel='x_only'):
    value = dict(channel=channel, steps=5, Ly=5, Lz=5, p=.02, shots=100,
                 target_shots=100, seed=123, point_seed='123456789abcdef0',
                 complete='true', elapsed_seconds=1.)
    value.update({f'joint_{i:02d}': n for i,n in enumerate((50,30,10,10))})
    for arm, count in [('final_only',40),('full_prefix',20)]:
        value[arm+'_failures']=count
        value[arm+'_logical_1_failures']=count
        value[arm+'_logical_2_failures']=0
    return value


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(M, 'Full-prefix analysis is not implemented')

    def test_paired_direction_and_intervals(self):
        result=M.paired([50,30,10,10])
        self.assertEqual(result['difference'],-.2)
        self.assertEqual(result['second'],'full_prefix')
        self.assertEqual(result['first_only'],30)
        self.assertEqual(result['second_only'],10)
        self.assertLess(result['ci_low'],-.2)
        self.assertGreater(result['ci_high'],-.2)

    def test_paired_endpoints_remain_bounded(self):
        for mask in range(4):
            bins=[0]*4; bins[mask]=100
            result=M.paired(bins)
            self.assertEqual(result['difference'], int(bool(mask&2))-int(bool(mask&1)))
            self.assertLessEqual(result['ci_low'],result['difference'])
            self.assertGreaterEqual(result['ci_high'],result['difference'])
            self.assertGreaterEqual(result['ci_low'],-1)
            self.assertLessEqual(result['ci_high'],1)

    def test_joint_logical_marginals_and_channel(self):
        for channel in ('x_only','z_only'):
            parsed=M.validate_point(row(channel))
            self.assertEqual(parsed['joint'],[50,30,10,10])
            self.assertEqual(parsed['channel'],channel)
        with self.assertRaises(ValueError): M.validate_point(row('mixed'))

    def test_corrupt_counts_are_rejected(self):
        for key,value in [('shots',99),('full_prefix_failures',21),
                          ('full_prefix_logical_1_failures',0),('joint_00',49),
                          ('complete','false'),('point_seed','1'),('Ly',7),
                          ('p',float('nan'))]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                M.validate_point(dict(row(),**{key:value}))

    def test_invalid_histogram_rejected(self):
        for bins in ([0]*4,[-1,1,1,1],[1,2,3],[1.2,2,3,4],[True,1,2,3]):
            with self.subTest(bins=bins), self.assertRaises(ValueError): M.paired(bins)

    def test_full_prefix_adapter_preserves_guarded_numerical_result(self):
        points=[]
        for size,counts in ((3,(2,4,6)),(5,(1,4,8))):
            for rate,count in zip((.01,.02,.03),counts):
                points.append(dict(steps=size,Ly=size,Lz=size,p=rate,shots=100,
                                   final_only_failures=30,full_prefix_failures=count))
        original=copy.deepcopy(points)
        legacy=[dict(steps=3,Ly=3,Lz=3,p=.01,shots=100,one_step_failures=2),
                dict(steps=3,Ly=3,Lz=3,p=.02,shots=100,one_step_failures=4),
                dict(steps=3,Ly=3,Lz=3,p=.03,shots=100,one_step_failures=6),
                dict(steps=5,Ly=5,Lz=5,p=.01,shots=100,one_step_failures=1),
                dict(steps=5,Ly=5,Lz=5,p=.02,shots=100,one_step_failures=4),
                dict(steps=5,Ly=5,Lz=5,p=.03,shots=100,one_step_failures=8)]
        actual=M.scaling(points,[.01,.02,.03],0,22)
        expected=M.OLD.scaling_diagnostics(legacy,[.01,.02,.03],0,22)
        self.assertEqual(actual,expected)
        self.assertEqual(points,original)
        self.assertEqual(actual['status'],'diagnostic_only')
        self.assertIsNone(actual['pc'])
        self.assertIsNone(actual['nu'])
        self.assertIn('grid_resolution',actual['failed_guards'])


class FigureTests(unittest.TestCase):
    def test_drop_smallest_diagnostic_is_saved_when_candidate_is_unavailable(self):
        spec=importlib.util.spec_from_file_location('full_prefix_drop_plot',PATH.with_name('plot_drop_smallest.py'))
        plot=importlib.util.module_from_spec(spec); spec.loader.exec_module(plot)
        report=dict(fits={'x_only':dict(sizes=[3,5],rates=[.01,.02],
            sensitivity=[dict(name='drop_smallest',candidate=None)])},marginals=[
            dict(channel='x_only',arm='full_prefix',L=5,p=.01,rate=0.,ci_low=0.,ci_high=.04),
            dict(channel='x_only',arm='full_prefix',L=5,p=.02,rate=1.,ci_low=.96,ci_high=1.)])
        with tempfile.TemporaryDirectory() as tmp:
            plot.plot(report,Path(tmp))
            for extension in ('png','pdf','svg'):
                path=Path(tmp)/f'x-full-prefix-drop-L3-fit.{extension}'
                self.assertTrue(path.exists(),f'Missing diagnostic: {path.name}')
                self.assertGreater(path.stat().st_size,1000)


def fixture_toml(data):
    def value(item):
        if isinstance(item,dict):
            return '{'+', '.join(json.dumps(k)+' = '+value(v) for k,v in item.items())+'}'
        if isinstance(item,list): return '['+', '.join(value(v) for v in item)+']'
        return json.dumps(item)
    return ''.join(k+' = '+value(v)+'\n' for k,v in data.items())


def fixture_protocol(channel):
    """Independent literal metadata from the approved runner contract."""
    x=channel=='x_only'
    return dict(circuit='legacy H+6-CNOT plaquette transfer',channel=channel,
        geometry='T=L; Ly=Lz=L',arm_order=['final_only','full_prefix'],
        initial_state='ideal encoded logical |00>' if x else 'noisy seed-H preparation of logical |++>',
        seed_preparation='ideal; zero seed fault locations' if x else '3 L^2 noisy seed-H outputs; no separate r=0 seed correction in either arm',
        final_only='no intermediate physical correction; exactly one final recovery and logical outcome',
        full_prefix='complete physical recovery on all available data slices Q0..Qr after EVERY transfer 1:T, including the last; complete-prefix syndrome zero; exactly one final recovery and logical outcome afterward',
        feedback_timing='after transfers 1:T inclusive; no separate r=0 correction',
        fault_locations='13 L^3 participating transfer outputs' if x else '3 L^2 seed-H outputs plus 13 L^3 participating transfer outputs',
        noise='independent Bernoulli-p faults in the selected channel at its recorded circuit locations',
        sampler='sample_toric3d_faults(p_x=p,p_z=0)' if x else 'sample_toric3d_diagonal_z_faults: seed first, then participating transfer outputs in schedule order',
        shared_final_decoder='family.prefixes[end]; full-plaquette BP+OSD with declared p prior; syndrome-validated recovery' if x else 'family.prefixes[end]; unit-weight data-edge PyMatching MWPM; syndrome-validated recovery',
        decoder_settings=dict(bp_method='minimum_sum',max_iter=50,ms_scaling_factor=1.,omp_thread_count=1,
                              osd_method='OSD_CS',osd_order=2,schedule='parallel') if x else dict(weights=1.,merge_strategy='smallest-weight'),
        z_weights='unit data-edge weights',failure='either residual logical parity after final recovery',
        joint_histogram_encoding='mask = final_only + 2*full_prefix; each term is its state-failure Boolean',
        point_order='sizes outer loop, rates inner loop; both strictly increasing',
        point_seed_generation='MersenneTwister(integer base seed); repeated scalar rand(rng,UInt64); discard duplicate values; 16 lowercase hexadecimal digits',
        point_rng='Random.Xoshiro(point_seed)',
        resume='recreate point Xoshiro and sample/discard completed records without decoding; rebuild raw CSV from authoritative checkpoint')


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(M, 'Full-prefix analysis is not implemented')

    def fixture(self,parent,name,channel,seed,start=100.):
        directory=Path(parent)/name; directory.mkdir()
        config=dict(channel=channel,sizes=[5],rates=[.01,.02],shots=100,seed=seed,batch_size=25,dimensions=[[5,5,5]])
        runtime=dict(julia_version='1.12.6',julia_threads=1,blas_threads=1,python_executable=sys.executable,
            python_version=sys.version,packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','ldpc','pymatching')},
            environment={key:'1' for key in ('JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
        runtime['environment'].update(JULIA_CONDAPKG_BACKEND='Null',JULIA_PYTHONCALL_EXE=sys.executable)
        sources=[M.ROOT/'Project.toml',PATH.with_name('run_study.jl'),*sorted((M.ROOT/'src').rglob('*.jl'))]
        if (M.ROOT/'Manifest.toml').exists(): sources.append(M.ROOT/'Manifest.toml')
        hashes={str(p.relative_to(M.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        rows=[dict(row(channel),p=p,seed=seed,point_seed=f'{seed+i:016x}') for i,p in enumerate(config['rates'])]
        states=[dict(completed_shots=100,elapsed_seconds=1.,joint_counts=[50,30,10,10],failures=[40,20],
                     logical_1_counts=[40,20],logical_2_counts=[0,0]) for _ in rows]
        checkpoint=dict(format_version=1,config=config,runtime=runtime,source_hashes=hashes,protocol=fixture_protocol(channel),
            point_seeds_hex=[r['point_seed'] for r in rows],states=states,generation=1,started_at_unix=start,updated_at_unix=start+1)
        self.seal(directory,rows,checkpoint)
        return directory,rows,checkpoint

    def seal(self,directory,rows,checkpoint):
        with (directory/'raw.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        checkpoint.pop('payload_sha256',None)
        checkpoint['payload_sha256']=hashlib.sha256(fixture_toml(checkpoint).encode()).hexdigest()
        blob=fixture_toml(checkpoint).encode(); (directory/'checkpoint.toml').write_bytes(blob)
        snapshot=f'checkpoints/checkpoint-{checkpoint["generation"]:06d}-{checkpoint["payload_sha256"][:12]}.toml'
        (directory/'checkpoints').mkdir(exist_ok=True); (directory/snapshot).write_bytes(blob)
        metadata={k:copy.deepcopy(v) for k,v in checkpoint.items() if k not in ('states','payload_sha256')}
        metadata.update(complete=True,completed_points=len(rows),expected_points=len(rows),completed_shots=sum(r['shots'] for r in rows),
                        raw_sha256=M.sha(directory/'raw.csv'),checkpoint_sha256=M.sha(directory/'checkpoint.toml'),current_snapshot=snapshot)
        (directory/'metadata.toml').write_text(fixture_toml(metadata))
        return metadata

    def production(self,parent):
        inputs=[]; checkpoints={}
        for channel,seed in (('x_only',100),('z_only',200)):
            directory,_,checkpoint=self.fixture(parent,'production-'+channel[0],channel,seed)
            inputs.append(directory); checkpoints[channel]=checkpoint
            self.fixture(parent,'pilot-'+channel[0],channel,seed+1000,start=1.)
        self.fixture(parent,'pilot-z-refine','z_only',1300,start=3.)
        points,records=M.load_points(inputs)
        numerical=('results/toric3d-h6-z-one-step-2026-10-15/analyze.py','examples/fit_toric3d_h6_threshold.py','examples/fit_toric3d_diagonal_feedback.py')
        selection=dict(stage='production',frozen_before_sampling=True,frozen_at_unix=10.,pilot_excluded=True,
            shards={'production-'+c[0]:{k:cp['config'][k] for k in ('channel','sizes','rates','shots','seed')} for c,cp in checkpoints.items()},
            channels={c:dict(sizes=[5],rates=[.01,.02],shots=100,fit_rates=[.01,.02],bootstrap_seed=50+i) for i,c in enumerate(checkpoints)},
            source_hashes=checkpoints['x_only']['source_hashes'],protocols={c:cp['protocol'] for c,cp in checkpoints.items()},
            numerical_source_sha256={n:M.sha(M.ROOT/n) for n in numerical},
            pilot_input_sha256={name:{f:M.sha(Path(parent)/name/f) for f in ('raw.csv','metadata.toml','checkpoint.toml')} for name in ('pilot-x','pilot-z')},
            additional_pilot_input_sha256={'pilot-z-refine':{f:M.sha(Path(parent)/'pilot-z-refine'/f) for f in ('raw.csv','metadata.toml','checkpoint.toml')}},
            bootstrap_replicates=200,total_fault_records=400)
        path=Path(parent)/'selection.json'; path.write_text(json.dumps(selection))
        return points,records,selection,path

    def test_complete_two_channel_provenance_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            points,records,selection,path=self.production(tmp)
            evidence=M.validate_selection(points,records,selection,path)
            self.assertTrue(evidence['validated'])
            self.assertEqual(evidence['additional_pilot_input_sha256'],selection['additional_pilot_input_sha256'])

    def test_refinement_pilot_manifest_is_required_and_complete(self):
        for mutation in ('absent','empty','null','missing_file','changed_hash','unexpected_pilot'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                points,records,selection,path=self.production(tmp)
                if mutation=='absent': selection.pop('additional_pilot_input_sha256')
                elif mutation=='empty': selection['additional_pilot_input_sha256']={}
                elif mutation=='null': selection['additional_pilot_input_sha256']=None
                elif mutation=='missing_file': selection['additional_pilot_input_sha256']['pilot-z-refine'].pop('checkpoint.toml')
                elif mutation=='changed_hash': selection['additional_pilot_input_sha256']['pilot-z-refine']['raw.csv']='0'*64
                else: selection['additional_pilot_input_sha256']['unfrozen-pilot']={}
                with self.assertRaises(ValueError): M.validate_selection(points,records,selection,path)

    def test_refinement_pilot_records_are_independently_validated(self):
        mutations=('after_freeze','wrong_channel','source','runtime','incomplete','histogram',
                   'production_seed_overlap','initial_pilot_seed_overlap')
        for mutation in mutations:
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                points,records,selection,path=self.production(tmp)
                directory=Path(tmp)/'pilot-z-refine'
                with (directory/'raw.csv').open(newline='') as stream: rows=list(csv.DictReader(stream))
                for row in rows: row['shots']=int(row['shots'])
                checkpoint=tomllib.loads((directory/'checkpoint.toml').read_text())
                if mutation=='after_freeze': checkpoint['updated_at_unix']=11.
                elif mutation=='wrong_channel':
                    checkpoint['config']['channel']='x_only'
                    checkpoint['protocol']=fixture_protocol('x_only')
                    for row in rows: row['channel']='x_only'
                elif mutation=='source': checkpoint['source_hashes']['Project.toml']='0'*64
                elif mutation=='runtime': checkpoint['runtime']['julia_version']='1.12.7'
                elif mutation=='incomplete': rows[0]['complete']='false'
                elif mutation=='histogram': checkpoint['states'][0]['joint_counts']=[49,31,10,10]
                else:
                    seed=100 if mutation=='production_seed_overlap' else 1100
                    rows[0]['point_seed']=f'{seed:016x}'
                    checkpoint['point_seeds_hex'][0]=rows[0]['point_seed']
                self.seal(directory,rows,checkpoint)
                selection['additional_pilot_input_sha256']['pilot-z-refine']={f:M.sha(directory/f) for f in ('raw.csv','metadata.toml','checkpoint.toml')}
                with self.assertRaises(ValueError): M.validate_selection(points,records,selection,path)

    def test_protocol_timing_seed_scope_and_decoder_changes_are_rejected(self):
        for channel in ('x_only','z_only'):
            for key,value in (('circuit','H+5'),('feedback_timing','before transfer'),('seed_preparation','separate seed repair'),
                              ('shared_final_decoder','different decoders'),('decoder_settings',{}),('full_prefix','skip final transfer')):
                with self.subTest(channel=channel,key=key),tempfile.TemporaryDirectory() as tmp:
                    directory,rows,checkpoint=self.fixture(tmp,'input',channel,10)
                    checkpoint['protocol'][key]=value; self.seal(directory,rows,checkpoint)
                    with self.assertRaises(ValueError): M.load_points([directory])

    def test_two_slice_provenance_cannot_be_relabelled_full_prefix(self):
        for channel in ('x_only','z_only'):
            with self.subTest(channel=channel),tempfile.TemporaryDirectory() as tmp:
                directory,rows,checkpoint=self.fixture(tmp,'input',channel,10)
                checkpoint['protocol']['arm_order']=['final_only','two_slice']
                checkpoint['protocol']['two_slice']=checkpoint['protocol'].pop('full_prefix')
                self.seal(directory,rows,checkpoint)
                with self.assertRaises(ValueError): M.load_points([directory])

    def test_source_raw_checkpoint_and_snapshot_integrity_are_enforced(self):
        for mutation in ('source','missing_source','raw_hash','checkpoint_histogram','snapshot','payload'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                directory,rows,checkpoint=self.fixture(tmp,'input','x_only',10)
                if mutation=='source': checkpoint['source_hashes']['Project.toml']='0'*64
                elif mutation=='missing_source': checkpoint['source_hashes'].pop('Project.toml')
                elif mutation=='checkpoint_histogram': checkpoint['states'][0]['joint_counts']=[49,31,10,10]
                metadata=self.seal(directory,rows,checkpoint)
                if mutation=='raw_hash': (directory/'raw.csv').write_text((directory/'raw.csv').read_text()+'\n')
                elif mutation=='snapshot': (directory/metadata['current_snapshot']).write_text('corrupt snapshot')
                elif mutation=='payload':
                    checkpoint['payload_sha256']='0'*64
                    (directory/'checkpoint.toml').write_text(fixture_toml(checkpoint))
                    metadata['checkpoint_sha256']=M.sha(directory/'checkpoint.toml')
                    (directory/'metadata.toml').write_text(fixture_toml(metadata))
                with self.assertRaises(ValueError): M.load_points([directory])

    def test_recorded_thread_backend_and_python_release_are_enforced(self):
        mutations=('OMP_NUM_THREADS','JULIA_CONDAPKG_BACKEND','JULIA_PYTHONCALL_EXE','python_version','julia_version')
        for key in mutations:
            with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
                directory,rows,checkpoint=self.fixture(tmp,'input','x_only',10)
                if key in ('python_version','julia_version'): checkpoint['runtime'][key]=''
                else: checkpoint['runtime']['environment'][key]='incorrect'
                self.seal(directory,rows,checkpoint)
                with self.assertRaises(ValueError): M.load_points([directory])

    def test_embedded_python_build_text_may_differ(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory,rows,checkpoint=self.fixture(tmp,'input','z_only',10)
            checkpoint['runtime']['python_version']=sys.version.split()[0]+' | separately compiled libpython'
            self.seal(directory,rows,checkpoint)
            self.assertEqual(len(M.load_points([directory])[0]),2)

    def test_missing_numerical_or_pilot_manifests_are_rejected(self):
        for mutation in ('numerical_empty','numerical_missing','pilot_empty','pilot_missing_shard','pilot_missing_file'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                points,records,selection,path=self.production(tmp)
                if mutation=='numerical_empty': selection['numerical_source_sha256']={}
                elif mutation=='numerical_missing': selection['numerical_source_sha256'].pop('examples/fit_toric3d_h6_threshold.py')
                elif mutation=='pilot_empty': selection['pilot_input_sha256']={}
                elif mutation=='pilot_missing_shard': selection['pilot_input_sha256'].pop('pilot-z')
                else: selection['pilot_input_sha256']['pilot-x'].pop('checkpoint.toml')
                with self.assertRaises(ValueError): M.validate_selection(points,records,selection,path)

    def test_frozen_channel_grid_and_bootstrap_settings_cannot_be_omitted(self):
        mutations=('channel','empty_fit','zero_fit','reverse_fit','duplicate_fit','bootstrap_zero','bootstrap_bool','negative_seed','boolean_seed','not_frozen')
        for mutation in mutations:
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as tmp:
                points,records,selection,path=self.production(tmp)
                if mutation=='channel': selection['channels'].pop('z_only')
                elif mutation=='empty_fit': selection['channels']['x_only']['fit_rates']=[]
                elif mutation=='zero_fit': selection['channels']['x_only']['fit_rates']=[0]
                elif mutation=='reverse_fit': selection['channels']['x_only']['fit_rates']=[.02,.01]
                elif mutation=='duplicate_fit': selection['channels']['x_only']['fit_rates']=[.01,.01]
                elif mutation=='bootstrap_zero': selection['bootstrap_replicates']=0
                elif mutation=='bootstrap_bool': selection['bootstrap_replicates']=True
                elif mutation=='negative_seed': selection['channels']['x_only']['bootstrap_seed']=-1
                elif mutation=='boolean_seed': selection['channels']['x_only']['bootstrap_seed']=True
                else: selection['frozen_before_sampling']=False
                with self.assertRaises(ValueError): M.validate_selection(points,records,selection,path)


if __name__=='__main__': unittest.main()
