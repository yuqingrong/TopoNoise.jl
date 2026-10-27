#!/usr/bin/env python3
"""Run frozen physical-x H+6 grouped shards with bounded local concurrency."""
import argparse
from collections import deque
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import tomllib

JULIA = '/Users/rongyuqing/.juliaup/bin/julia'
PYTHON = '/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python'
THREAD_KEYS = ('JULIA_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(text)
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_checkpoint(path):
    raw = Path(path).read_bytes()
    checkpoint = tomllib.loads(raw.decode())
    payload, count = re.subn(rb'^payload_sha256 = "[0-9a-f]{64}"\n', b'', raw, flags=re.M)
    if count != 1 or hashlib.sha256(payload).hexdigest() != checkpoint.get('payload_sha256'):
        raise ValueError(f'checkpoint checksum mismatch: {path}')
    return checkpoint


def mirrors_current(directory):
    """A committed final checkpoint may still need its CSV mirrors repaired."""
    directory = Path(directory)
    try:
        checkpoint_path = directory/'checkpoint.toml'
        checkpoint = read_checkpoint(checkpoint_path)
        metadata = tomllib.loads((directory/'metadata.toml').read_text())
        return (metadata.get('generation') == checkpoint['generation']
                and metadata.get('config') == checkpoint['config']
                and metadata.get('checkpoint_sha256') == digest(checkpoint_path)
                and metadata.get('raw_sha256') == digest(directory/'raw.csv'))
    except (OSError, ValueError, KeyError):
        return False


def signal_worker(process, sig):
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass  # The worker can exit between poll() and killpg().


def validate_grid(selection):
    if selection['shots'] != 1000:
        raise ValueError('production requires exactly 1000 shots per point')
    expected = {(ch, size, p) for ch, config in selection['channels'].items()
                for size in config['sizes'] for p in config['rates']}
    actual, seeds, point_seeds = [], [], []
    for name, shard in selection['shards'].items():
        if Path(name).name != name or name in ('.', '..'):
            raise ValueError('unsafe shard name')
        rates=shard['rates']
        if rates != sorted(set(rates)):
            raise ValueError('shard grid rates must be unique and increasing')
        actual.extend((shard['channel'], shard['size'], p) for p in rates)
        seeds.append(shard['seed'])
        stream_seeds=shard['point_seeds_hex']
        if len(stream_seeds) != len(rates) or any(not re.fullmatch(r'[0-9a-f]{16}', v) for v in stream_seeds):
            raise ValueError('point seed count or format mismatch')
        point_seeds.extend(stream_seeds)
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise ValueError('shard grid differs from frozen selection')
    if len(seeds) != len(set(seeds)) or len(point_seeds) != len(set(point_seeds)):
        raise ValueError('sampling seeds overlap')
    if len(actual) * selection['shots'] != selection['total_fault_records']:
        raise ValueError('sample total mismatch')
    if selection['total_protocol_outcomes'] != 3 * selection['total_fault_records']:
        raise ValueError('outcome total mismatch')


def verify_frozen(study, selection):
    validate_grid(selection)
    if selection.get('stage') != 'direction_comparison' or not selection.get('frozen_before_sampling'):
        raise ValueError('selection has not been frozen for the direction comparison')
    archive = (study/'archive').resolve()
    for relative, sha in {**selection['source_hashes'], **selection['auxiliary_source_hashes']}.items():
        path=(archive/relative).resolve()
        if not path.is_relative_to(archive) or digest(path) != sha:
            raise ValueError(f'frozen source changed: {relative}')
    for kind in ('csv','report'):
        if digest(selection['reference'][kind]) != selection['reference'][kind+'_sha256']:
            raise ValueError(f'diagonal reference changed: {kind}')
    for filename,key in (('verification.json','verification_sha256'),('timing-pilot.toml','timing_pilot_sha256')):
        if digest(archive/filename) != selection[key]:
            raise ValueError(f'frozen pre-production evidence changed: {filename}')


def ordered_jobs(selection):
    """Interleave X/Z initially; independent completed shards can be skipped."""
    by_channel = {ch: deque(sorted((name for name, item in selection['shards'].items()
                        if item['channel'] == ch),
                        key=lambda name: (selection['shards'][name]['size'],
                                          selection['shards'][name]['rates'][0])))
                  for ch in selection['channels']}
    jobs = []
    while any(by_channel.values()):
        for queue in by_channel.values():
            if queue:
                jobs.append(queue.popleft())
    return jobs


def runner_command(study, selection, name, resume=False):
    shard = selection['shards'][name]
    archive = Path(study)/'archive'
    command = [JULIA, '--compiled-modules=existing', f'--project={archive / "runtime"}',
               '--startup-file=no', str(archive/'run_study.jl'), '--channel', shard['channel'],
               '--output', str(archive/'runs'/name), '--sizes', str(shard['size']),
               '--rates', ','.join(map(repr, shard['rates'])), '--shots', str(selection['shots']),
               '--seed', str(shard['seed']), '--batch-size', str(selection['batch_size'])]
    if resume:
        command.append('--resume')
    return command


def environment():
    env = dict(os.environ)
    env.update({name: '1' for name in THREAD_KEYS})
    env.update(JULIA_CONDAPKG_BACKEND='Null', JULIA_PYTHONCALL_EXE=PYTHON,
               MPLCONFIGDIR='/private/tmp/toponoise-h6-mpl')
    return env


def collect_status(study, selection, active, phase, started, error=None):
    points, completed, completed_points = {}, 0, 0
    for name, shard in selection['shards'].items():
        path = study/'archive/runs'/name/'checkpoint.toml'
        counts = [0]*len(shard['rates'])
        if path.exists():
            checkpoint = read_checkpoint(path)
            config = checkpoint['config']
            expected = {'channel': shard['channel'], 'sizes': [shard['size']],
                        'rates': shard['rates'], 'shots': selection['shots'], 'seed': shard['seed']}
            if any(config[k] != v for k, v in expected.items()) or len(checkpoint['states']) != len(counts):
                raise ValueError(f'shard configuration mismatch: {name}')
            if checkpoint.get('point_seeds_hex') != shard['point_seeds_hex']:
                raise ValueError(f'shard point seeds mismatch: {name}')
            counts = [state['completed_shots'] for state in checkpoint['states']]
            if any(type(n) is not int or not 0 <= n <= selection['shots'] for n in counts):
                raise ValueError(f'invalid shot count: {name}')
        n=sum(counts)
        done=sum(k == selection['shots'] for k in counts)
        completed += n
        completed_points += done
        points[name] = {'completed_shots': n, 'complete': done == len(counts),
                        'completed_points': done, 'target_points':len(counts), 'counts':counts}
    return {'phase': phase, 'supervisor_pid': os.getpid(), 'started_at_unix': started,
            'updated_at_unix': time.time(), 'completed_shots': completed,
            'target_shots': selection['total_fault_records'],
            'completed_protocol_outcomes': 3*completed,
            'completed_points': completed_points,
            'target_points': sum(len(s['rates']) for s in selection['shards'].values()),
            'active': {name: p.pid for name, (p, log) in active.items()},
            'error': error, 'points': points}


def analyze(study, progress):
    command = [PYTHON, str(study/'archive/analyze.py'), '--study', str(study)]
    if progress:
        command.append('--progress')
    with (study/'archive/analysis.log').open('a') as log:
        subprocess.run(command, cwd=study, env=environment(), stdout=log,
                       stderr=subprocess.STDOUT, check=True)


def supervise(study, workers):
    study = Path(study).resolve()
    selection = json.loads((study/'data/selection.json').read_text())
    verify_frozen(study, selection)
    archive = study/'archive'
    lock = (archive/'supervisor.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError('A supervisor already owns this study')
    active, started, last_export = {}, time.time(), 0.0
    status_path = study/'data/progress.json'
    def save(phase, error=None):
        try:
            state = collect_status(study, selection, active, phase, started, error)
        except (OSError, ValueError, KeyError) as status_error:
            if phase not in ('failed', 'stopped'):
                raise
            # A corrupt checkpoint must not hide the original failure status.
            try:
                state = json.loads(status_path.read_text())
            except (OSError, ValueError):
                state = {}
            state.update(phase=phase, error=error, status_error=str(status_error),
                         updated_at_unix=time.time(), supervisor_pid=os.getpid(), active={})
        atomic_text(status_path, json.dumps(state, indent=2)+'\n')
        return state
    stopping = False
    def stop(signum, frame):
        nonlocal stopping
        stopping = True
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, stop)
    # Keeps an explicitly launched long simulation running during idle time;
    # closing the laptop or shutting it down still suspends/ends computation.
    keep_awake = None
    if Path('/usr/bin/caffeinate').exists():
        keep_awake = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())])
    try:
        initial = save('running')
        jobs = deque(name for name in ordered_jobs(selection)
                     if not initial['points'][name]['complete']
                     or not mirrors_current(archive/'runs'/name))
        while jobs or active:
            if stopping:
                raise InterruptedError('Stopped; committed batches can be resumed')
            changed = False
            for name, (process, log) in list(active.items()):
                code = process.poll()
                if code is None:
                    continue
                log.close()
                del active[name]
                if code != 0:
                    raise RuntimeError(f'{name} exited with code {code}; inspect archive/runs/{name}/run.log')
                check = collect_status(study, selection, active, 'running', started)
                if not check['points'][name]['complete']:
                    raise RuntimeError(f'{name} exited without completing its target')
                changed = True
                print(f'Completed {name}', flush=True)
            while jobs and len(active) < workers:
                name = jobs.popleft()
                directory = archive/'runs'/name
                directory.mkdir(parents=True, exist_ok=True)
                resume = (directory/'checkpoint.toml').exists()
                log = (directory/'run.log').open('a')
                try:
                    process = subprocess.Popen(runner_command(study, selection, name, resume),
                        cwd=archive/'runtime', env=environment(), stdout=log,
                        stderr=subprocess.STDOUT, start_new_session=True,
                        pass_fds=(lock.fileno(),))
                except BaseException:
                    log.close()
                    raise
                active[name] = process, log
                print(f'{"Resumed" if resume else "Started"} {name} pid={process.pid}', flush=True)
            state = save('running')
            # Consolidated exports are refreshed at most once per five minutes.
            if 0 < state['completed_shots'] < state['target_shots'] and time.time()-last_export >= 300:
                analyze(study, progress=True)
                last_export = time.time()
            if active:
                time.sleep(10)
        save('reporting')
        analyze(study, progress=False)
        save('complete')
        print('All physical-axis samples and the direction-comparison PDF are complete.', flush=True)
    except BaseException as error:
        for process, log in active.values():
            if process.poll() is None:
                signal_worker(process, signal.SIGTERM)
        for process, log in active.values():
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                signal_worker(process, signal.SIGKILL)
                process.wait()
            log.close()
        active.clear()
        save('stopped' if isinstance(error, InterruptedError) else 'failed', str(error))
        raise
    finally:
        if keep_awake is not None:
            keep_awake.terminate()
        lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--workers', type=int, default=2, choices=(1, 2, 3, 4))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    selection = json.loads((args.study/'data/selection.json').read_text())
    if args.dry_run:
        verify_frozen(args.study, selection)
        print(json.dumps({'points': sum(len(s['rates']) for s in selection['shards'].values()), 'shots': selection['shots'],
                          'workers': args.workers, 'records': selection['total_fault_records']}))
    else:
        supervise(args.study, args.workers)


if __name__ == '__main__':
    main()
