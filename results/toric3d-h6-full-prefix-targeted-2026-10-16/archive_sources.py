#!/usr/bin/env python3
"""Plan, smoke-check, and (only after final audits) archive this study's sources.

No production data, fitted results, figures, environments, or prior data bundles
are included. Member names preserve paths relative to the project root.
"""
import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REL = HERE.relative_to(ROOT)
PREVIOUS = Path('results/toric3d-h6-full-prefix-reanalysis-2026-10-16')
OLD = Path('results/toric3d-h6-full-prefix-2026-10-15')
MANIFEST = 'SOURCE_MANIFEST.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_relative(value):
    path = Path(value)
    require(not path.is_absolute() and '..' not in path.parts and str(path) not in ('', '.'),
            f'Unsafe member path: {value}')
    return path


def source_plan():
    plan = json.loads((HERE / 'selection.json').read_text())
    require(sha(HERE / 'selection.json') == (HERE / 'selection.sha256').read_text().split()[0],
            'Frozen selection digest differs')
    frozen = plan['simulation_source_sha256']
    require(len(frozen) == 56, 'Expected 56 frozen simulation source files')
    for name, digest in frozen.items():
        require(sha(ROOT / checked_relative(name)) == digest, f'Frozen source changed: {name}')
    setup = json.loads((HERE / 'setup-freeze.json').read_text())
    for name, digest in setup['files_sha256'].items():
        require(sha(HERE / checked_relative(name)) == digest, f'Frozen setup changed: {name}')
    prior = {}
    for name, digest in plan['prior_analysis_report_sha256'].items():
        path = ROOT / PREVIOUS / name
        require(sha(path) == digest, f'Prior report changed: {name}')
        prior[name] = json.loads(path.read_text())
    numerical = prior['x-report.json']['provenance']['numerical_source_sha256']
    require(set(numerical) == {'examples/fit_toric3d_h6_threshold.py',
                               'examples/fit_toric3d_diagonal_feedback.py'},
            'Transitive numerical dependency list changed; review the archive plan')
    numerical = {**numerical,
        str(PREVIOUS / 'x_analysis.py'): prior['x-report.json']['provenance']['analysis_sha256'],
        str(PREVIOUS / 'z_models.py'): prior['z-report.json']['source_sha256']}
    for name, digest in numerical.items():
        require(sha(ROOT / name) == digest, f'Numerical dependency changed: {name}')
    require(sha(ROOT / OLD / 'ARTIFACTS.sha256.json') == plan['old_manifest_sha256'],
            'Prior-study artifact manifest differs')
    old_manifest = json.loads((ROOT / OLD / 'ARTIFACTS.sha256.json').read_text())
    replay = 'verification/replay_point.jl'
    require(sha(ROOT / OLD / replay) == old_manifest[replay], 'Independent replay source changed')
    preserved = json.loads((ROOT / OLD / 'verification/prior-study-integrity.json').read_text())
    require(sha(ROOT / 'document/note.tex') == preserved['note_sha256'], 'Preserved note changed')

    groups = {
        'frozen_simulation': list(frozen),
        'analysis_dependencies': list(numerical),
        'plot_theme': ['examples/plotting_theme.py', 'examples/styles/paper-no-grid.mplstyle'],
        'independent_replay': [str(OLD / replay)],
        'project_and_preservation': ['CondaPkg.toml', 'LICENSE', 'document/note.tex',
                                     'test/toric3d/h6_full_prefix.jl'],
        'targeted_study': [str(REL / name) for name in (
            'selection.json', 'selection.sha256', 'jobs.toml', 'setup-freeze.json',
            'run_worker.jl', 'analyze.py', 'test_analysis.py', 'plot_targeted.py',
            'archive_sources.py', 'reproduce.txt', 'README.md', 'verification/preflight.toml')],
        'targeted_verification': [str(path.relative_to(ROOT))
            for pattern in ('*.py', '*.jl') for path in sorted((HERE / 'verification').glob(pattern))],
    }
    paths = sorted({name for group in groups.values() for name in group})
    files = {}
    for name in paths:
        path = ROOT / checked_relative(name)
        require(path.is_file() and not path.is_symlink(), f'Expected regular source file: {name}')
        files[name] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    return {
        'schema_version': 1,
        'archive_kind': 'source-only; project-relative paths; external data prerequisites required',
        'study_path': str(REL),
        'selection_sha256': sha(HERE / 'selection.json'),
        'groups': groups, 'file_count': len(files), 'files': files,
        'external_prerequisites': {
            'old_study': {'path': str(OLD), 'manifest': 'ARTIFACTS.sha256.json',
                'manifest_sha256': plan['old_manifest_sha256'],
                'required_scope': 'Complete separately retained old-study bundle, including every '
                                  'manifest-listed file; the preservation audit checks them all.'},
            'prior_reanalysis_reports': {str(PREVIOUS / name): digest
                                        for name, digest in plan['prior_analysis_report_sha256'].items()},
            'targeted_data': {'path': str(REL),
                'required_scope': 'Separately retained benchmark-x/ and benchmark-z/ fixtures; '
                    'all 13 production shards with raw.csv, metadata.toml, checkpoint.toml and '
                    'every immutable checkpoints/ snapshot; completed analysis outputs; replay fixtures.'},
            'runtime': 'Julia 1.12.6; Python 3.12.13; numpy 2.5.1; scipy 1.18.0; '
                       'ldpc 2.4.1; pymatching 2.4.0; matplotlib 3.11.1. Exact resume/replay '
                       'also requires the recorded Python executable path and one-thread environment.',
        },
        'exclusions': ['production and benchmark data', 'analysis outputs', 'figures',
                       'checkpoint archives', 'prior-study data bundles', 'installed runtimes',
                       '__pycache__', '.git', 'logs'],
    }


SMOKE_PROGRAM = r'''
import importlib.util, json, pathlib, sys
root=pathlib.Path(sys.argv[1]).resolve(); relative=pathlib.Path(sys.argv[2])
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
a=load('extracted_analysis',root/relative/'analyze.py')
p=load('extracted_plot',root/relative/'plot_targeted.py')
assert a.ROOT==p.ROOT==root
locations=[a.__file__,a.X.__file__,a.Z.__file__,a.X.H6.__file__,a.X.H6.ENGINE.__file__,
           p.__file__,sys.modules['plotting_theme'].__file__]
assert all(pathlib.Path(v).resolve().is_relative_to(root) for v in locations)
assert (root/'examples/styles/paper-no-grid.mplstyle').is_file()
for filename in ['check_production.py','check_numerics.py','test_check_production.py',
                 'test_check_numerics.py','test_archive_sources.py']:
    load('extracted_'+filename.removesuffix('.py'),root/relative/'verification'/filename)
load('extracted_analysis_tests',root/relative/'test_analysis.py')
plan=a.read_plan(root/relative)
assert len(a.fit_schedule(plan,'x_only'))==24 and len(a.fit_schedule(plan,'z_only'))==36
assert len(a.profile_schedule(plan))==12
lo,hi=a.X.clopper_pearson(5,10,.05);assert 0<lo<.5<hi<1
assert p.wilson_interval(0,100)[0]==0 and p.wilson_interval(100,100)[1]==1
with p.paper_theme():
    fig,ax=p.plt.subplots();p.style_axes(ax)
    assert not any(g.get_visible() for g in ax.get_xgridlines()+ax.get_ygridlines())
    p.plt.close(fig)
print(json.dumps({'status':'passed','imported_source_paths':[str(pathlib.Path(v).relative_to(root)) for v in locations],
                  'scope':'Isolated imports, schedule and arithmetic checks, and no-grid theme; '
                          'no old data, circuit simulation, fitting, production data, or PDF output required.'}))
'''


def smoke_tree(root, manifest):
    for name, item in manifest['files'].items():
        require(sha(root / checked_relative(name)) == item['sha256'], f'Extracted source differs: {name}')
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                       MPLCONFIGDIR=os.environ.get('MPLCONFIGDIR', str(root / '.smoke-matplotlib')),
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    result = subprocess.run([sys.executable, '-I', '-B', '-c', SMOKE_PROGRAM,
                             str(root), manifest['study_path']], cwd=root, env=environment,
                            check=True, text=True, capture_output=True, timeout=120)
    return json.loads(result.stdout)


def smoke_sources(manifest):
    with tempfile.TemporaryDirectory(prefix='h6-targeted-source-smoke-') as directory:
        root = Path(directory)
        for name in manifest['files']:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, path)
        return smoke_tree(root, manifest)


def verify_archive(path):
    with tempfile.TemporaryDirectory(prefix='h6-targeted-extract-smoke-') as directory:
        root = Path(directory)
        with tarfile.open(path, 'r:gz') as archive:
            members = archive.getmembers()
            require(len({m.name for m in members}) == len(members), 'Duplicate archive members')
            for member in members:
                require(member.isfile(), 'Only regular files are permitted in the source archive')
                destination = root / checked_relative(member.name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as source:
                    destination.write_bytes(source.read())
        manifest = json.loads((root / MANIFEST).read_text())
        require({m.name for m in members} == set(manifest['files']) | {MANIFEST},
                'Archive membership differs from its manifest')
        result = smoke_tree(root, manifest)
        return {'archive_sha256': sha(path), 'file_count': manifest['file_count'], 'smoke': result}


def require_final_audits():
    report_path = HERE / 'analysis/report.json'
    require(report_path.is_file(), 'Final archive creation requires a completed analysis report')
    report = json.loads(report_path.read_text())
    require(report.get('status') in ('complete', 'complete_with_numerical_failures')
            and report.get('completed_at_unix'), 'Analysis is not complete')
    for name, digest in report['provenance']['numerical_sources_sha256'].items():
        require(sha(ROOT / checked_relative(name)) == digest,
                f'Numerical source changed after the final analysis: {name}')
    audit_hashes = {}
    for name in ('production-audit.json', 'numerical-audit.json'):
        path = HERE / 'verification' / name
        require(path.is_file(), f'Final archive creation requires verification/{name}')
        audit = json.loads(path.read_text())
        require(audit.get('status') == 'passed', f'Final audit has not passed: {name}')
        require(audit.get('selection_sha256') == sha(HERE / 'selection.json'),
                f'Audit selection binding differs: {name}')
        audit_hashes[str(path.relative_to(ROOT))] = sha(path)
        if name == 'production-audit.json':
            plan = json.loads((HERE / 'selection.json').read_text())
            wanted = {job['name'] for job in plan['jobs']}
            shards = audit.get('shards', [])
            require(len(shards) == len(wanted) and {shard['name'] for shard in shards} == wanted,
                    'Production audit does not bind every planned shard exactly once')
            for shard in shards:
                directory = HERE / checked_relative(shard['name'])
                for filename, key in (('raw.csv', 'raw_sha256'),
                                      ('checkpoint.toml', 'checkpoint_sha256'),
                                      ('metadata.toml', 'metadata_sha256')):
                    require(shard.get(key) == sha(directory / filename),
                            f'Production audit is stale: {shard["name"]}/{filename}')
            if 'checker_sha256' in audit:
                require(audit['checker_sha256'] == sha(HERE / 'verification/check_production.py'),
                        'Production audit checker changed after the final audit')
        if name == 'numerical-audit.json':
            require(audit.get('report_sha256') == sha(report_path), 'Numerical audit is stale')
            require(audit.get('input_sha256') == sha(HERE / 'analysis/input.json')
                    == report['provenance'].get('input_sha256'), 'Numerical input binding is stale')
            require(audit.get('checker_sha256') == sha(HERE / 'verification/check_numerics.py'),
                    'Numerical audit checker changed after the final audit')
    return {'final_report_sha256': sha(report_path), 'final_audits_sha256': audit_hashes}


def create_archive(path, manifest):
    final_bindings = require_final_audits()
    manifest = {**manifest, **final_bindings}
    require(not path.exists(), f'Refusing to replace existing archive: {path}')
    receipt_path = path.with_name(path.name + '.receipt.json')
    require(not receipt_path.exists(), f'Refusing to replace existing receipt: {receipt_path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    payloads = {name: (ROOT / name).read_bytes() for name in manifest['files']}
    for name, payload in payloads.items():
        require(hashlib.sha256(payload).hexdigest() == manifest['files'][name]['sha256'],
                f'Source changed during packaging: {name}')
    payloads[MANIFEST] = (json.dumps(manifest, indent=2, sort_keys=True) + '\n').encode()
    temporary = path.with_name(path.name + '.tmp')
    require(not temporary.exists(), f'Temporary archive already exists: {temporary}')
    try:
        with temporary.open('xb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', filename='', mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w') as archive:
                for name, payload in sorted(payloads.items()):
                    member = tarfile.TarInfo(name)
                    member.size = len(payload); member.mode = 0o644; member.mtime = 0
                    archive.addfile(member, io.BytesIO(payload))
        checked = verify_archive(temporary)
        require(not path.exists(), f'Archive appeared during packaging: {path}')
        temporary.rename(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    receipt = {'archive': str(path.resolve()), **final_bindings, **checked}
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-file', type=Path, help='Write a provisional source manifest; no archive')
    parser.add_argument('--smoke', action='store_true', help='Copy only sources to a temporary tree and import them')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--create', type=Path, help='Create a source archive after completion and both final audits')
    modes.add_argument('--verify-archive', type=Path, help='Extract, hash-check and import an existing source archive')
    args = parser.parse_args()
    if args.verify_archive:
        print(json.dumps(verify_archive(args.verify_archive), indent=2)); return
    manifest = source_plan()
    if args.plan_file:
        args.plan_file.parent.mkdir(parents=True, exist_ok=True)
        args.plan_file.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    result = {'status': 'planned', 'file_count': manifest['file_count'],
              'frozen_simulation_files': len(manifest['groups']['frozen_simulation']),
              'source_bytes': sum(item['bytes'] for item in manifest['files'].values()),
              'archive_created': False}
    if args.smoke:
        result['smoke'] = smoke_sources(manifest)
    if args.create:
        result = {'status': 'created', 'archive_created': True, **create_archive(args.create, manifest)}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
