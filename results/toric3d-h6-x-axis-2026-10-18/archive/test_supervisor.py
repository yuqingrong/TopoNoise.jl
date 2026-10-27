import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

PATH = Path(__file__).with_name('supervise.py')
SPEC = importlib.util.spec_from_file_location('matched_supervisor', PATH)
M = importlib.util.module_from_spec(SPEC)
if PATH.exists():
    SPEC.loader.exec_module(M)


class SupervisorTests(unittest.TestCase):
    def setUp(self):
        self.selection = {
            'stage': 'direction_comparison', 'frozen_before_sampling': True,
            'shots':1000, 'batch_size':25, 'total_fault_records':4000,
            'total_protocol_outcomes':12000,
            'channels': {ch:{'sizes':[3], 'rates':[.01,.02]} for ch in ('x_only','z_only')},
            'shards': {ch+'-L3':{'channel':ch,'size':3,'rates':[.01,.02],
                'seed':i+10,'point_seeds_hex':[f'{i*2+j:016x}' for j in (1,2)]}
                for i,ch in enumerate(('x_only','z_only'))}}
        self.assertTrue(PATH.exists(), 'grouped physical-axis supervisor is not implemented')

    def test_jobs_cover_grid_and_start_both_channels(self):
        jobs = M.ordered_jobs(self.selection)
        self.assertEqual(len(jobs), len(self.selection['shards']))
        self.assertEqual(set(jobs), set(self.selection['shards']))
        self.assertEqual({self.selection['shards'][j]['channel'] for j in jobs[:2]}, {'x_only', 'z_only'})

    def test_command_preserves_exact_sampling_config(self):
        name = M.ordered_jobs(self.selection)[0]
        shard = self.selection['shards'][name]
        command = M.runner_command(PATH.parents[1], self.selection, name, resume=True)
        for flag, value in [('--shots', '1000'), ('--sizes', str(shard['size'])),
                            ('--rates', ','.join(map(repr,shard['rates']))), ('--seed', str(shard['seed']))]:
            self.assertEqual(command[command.index(flag)+1], value)
        self.assertIn('--resume', command)
        self.assertIn('--compiled-modules=existing', command)

    def test_atomic_write_and_checksum_reject_corruption(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'checkpoint.toml'
            payload = 'generation = 1\n\n[[states]]\ncompleted_shots = 5\n'
            digest = hashlib.sha256(payload.encode()).hexdigest()
            text = f'payload_sha256 = "{digest}"\n' + payload
            M.atomic_text(path, text)
            self.assertEqual(M.read_checkpoint(path)['states'][0]['completed_shots'], 5)
            path.write_text(text.replace('completed_shots = 5', 'completed_shots = 6'))
            with self.assertRaisesRegex(ValueError, 'checksum'):
                M.read_checkpoint(path)

    def test_selection_mismatch_refused(self):
        broken = copy.deepcopy(self.selection)
        first = next(iter(broken['shards']))
        broken['shards'][first]['rates'][0] = .49
        with self.assertRaisesRegex(ValueError, 'grid'):
            M.validate_grid(broken)

    def test_completed_checkpoint_needs_mirror_repair_after_interrupted_write(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            payload = 'generation = 3\n\n[config]\nshots = 1000\n\n[[states]]\ncompleted_shots = 1000\n'
            sha = hashlib.sha256(payload.encode()).hexdigest()
            (root/'checkpoint.toml').write_text(f'payload_sha256 = "{sha}"\n' + payload)
            self.assertFalse(M.mirrors_current(root))
            (root/'raw.csv').write_text('shots\n1000\n')
            (root/'metadata.toml').write_text(
                f'generation = 3\ncheckpoint_sha256 = "{M.digest(root/"checkpoint.toml")}"\n'
                f'raw_sha256 = "{M.digest(root/"raw.csv")}"\n\n[config]\nshots = 1000\n')
            self.assertTrue(M.mirrors_current(root))
            (root/'raw.csv').write_text('shots\n975\n')
            self.assertFalse(M.mirrors_current(root))

    def test_cleanup_tolerates_worker_exit_race(self):
        with mock.patch.object(M.os, 'killpg', side_effect=ProcessLookupError):
            M.signal_worker(mock.Mock(pid=123), M.signal.SIGTERM)

    def test_worker_keeps_lock_after_supervisor_closes_its_copy(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'lock'
            owner = path.open('a')
            M.fcntl.flock(owner, M.fcntl.LOCK_EX | M.fcntl.LOCK_NB)
            child = M.subprocess.Popen([M.sys.executable, '-c', 'import time; time.sleep(10)'],
                                       pass_fds=(owner.fileno(),))
            owner.close()
            try:
                with path.open('a') as contender:
                    with self.assertRaises(BlockingIOError):
                        M.fcntl.flock(contender, M.fcntl.LOCK_EX | M.fcntl.LOCK_NB)
            finally:
                child.terminate()
                child.wait()
            with path.open('a') as contender:
                M.fcntl.flock(contender, M.fcntl.LOCK_EX | M.fcntl.LOCK_NB)

    def test_corruption_still_writes_failed_status(self):
        with tempfile.TemporaryDirectory() as folder:
            study = Path(folder)
            (study/'data').mkdir()
            (study/'archive').mkdir()
            (study/'data/selection.json').write_text(json.dumps(self.selection))
            with mock.patch.object(M, 'verify_frozen'), \
                 mock.patch.object(M, 'collect_status', side_effect=ValueError('checksum mismatch')), \
                 mock.patch.object(M.subprocess, 'Popen'):
                with self.assertRaisesRegex(ValueError, 'checksum'):
                    M.supervise(study, 2)
            status = json.loads((study/'data/progress.json').read_text())
            self.assertEqual(status['phase'], 'failed')
            self.assertIn('checksum', status['error'])
            self.assertEqual(status['active'], {})
        broken = copy.deepcopy(self.selection)
        broken['shots'] = 9999
        with self.assertRaisesRegex(ValueError, '1000'):
            M.validate_grid(broken)


    def test_grouped_shards_cover_each_point_exactly_once(self):
        M.validate_grid(self.selection)
        for mutator in (
            lambda s: s['shards']['x_only-L3']['rates'].append(.01),
            lambda s: s['shards']['x_only-L3']['point_seeds_hex'].pop(),
            lambda s: s['shards']['x_only-L3']['point_seeds_hex'].__setitem__(0, s['shards']['z_only-L3']['point_seeds_hex'][0])):
            broken=copy.deepcopy(self.selection); mutator(broken)
            with self.assertRaises(ValueError): M.validate_grid(broken)

    def test_verified_evidence_cannot_change_after_freeze(self):
        with tempfile.TemporaryDirectory() as folder:
            study=Path(folder); archive=study/'archive'; archive.mkdir()
            for name in ('source.jl','analyze.py','verification.json','timing-pilot.toml','reference.csv','reference.json'):
                (archive/name).write_text(name)
            selection=copy.deepcopy(self.selection)
            selection.update(source_hashes={'source.jl':M.digest(archive/'source.jl')},
                auxiliary_source_hashes={'analyze.py':M.digest(archive/'analyze.py')},
                reference={kind:str(archive/name) for kind,name in (('csv','reference.csv'),('report','reference.json'))},
                verification_sha256=M.digest(archive/'verification.json'),
                timing_pilot_sha256=M.digest(archive/'timing-pilot.toml'))
            for kind in ('csv','report'):
                selection['reference'][kind+'_sha256']=M.digest(selection['reference'][kind])
            M.verify_frozen(study,selection)
            for name in ('verification.json','timing-pilot.toml'):
                path=archive/name; original=path.read_text(); path.write_text('changed')
                with self.assertRaisesRegex(ValueError,'evidence changed'):
                    M.verify_frozen(study,selection)
                path.write_text(original)

    def test_grouped_status_counts_points_not_only_jobs(self):
        with tempfile.TemporaryDirectory() as folder:
            study=Path(folder); run=study/'archive/runs/x_only-L3'; run.mkdir(parents=True)
            payload = ('generation = 1\npoint_seeds_hex = ["0000000000000001", "0000000000000002"]\n\n[config]\nchannel = "x_only"\nsizes = [3]\n'
                'rates = [0.01, 0.02]\nshots = 1000\nseed = 10\n\n'
                '[[states]]\ncompleted_shots = 1000\n\n[[states]]\ncompleted_shots = 25\n')
            sha=hashlib.sha256(payload.encode()).hexdigest()
            (run/'checkpoint.toml').write_text(f'payload_sha256 = "{sha}"\n'+payload)
            status=M.collect_status(study,self.selection,{},'running',1)
            self.assertEqual(status['completed_points'],1)
            self.assertEqual(status['target_points'],4)
            self.assertEqual(status['completed_shots'],1025)
            self.assertEqual(status['target_shots'],4000)
            self.assertFalse(status['points']['x_only-L3']['complete'])

if __name__ == '__main__':
    unittest.main()
