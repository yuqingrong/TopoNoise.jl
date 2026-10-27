"""Record the pilot-informed production design before any production sampling."""
from pathlib import Path
import hashlib
import json
import time
import tomllib
import tarfile

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
assert not (STUDY / 'selection.json').exists(), 'The frozen selection must not be overwritten.'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def input_hashes(folder):
    return {name: sha(STUDY / folder / name)
            for name in ('raw.csv', 'metadata.toml', 'checkpoint.toml')}


sources = json.loads((STUDY / 'verification/simulation-source-sha256.json').read_text())
assert all(sha(ROOT / path) == digest for path, digest in sources.items())
metas = {channel: tomllib.loads((STUDY / f'pilot-{channel[0]}/metadata.toml').read_text())
         for channel in ('x_only', 'z_only')}
refine = tomllib.loads((STUDY / 'pilot-z-refine/metadata.toml').read_text())
assert all(m['complete'] and m['source_hashes'] == sources for m in [*metas.values(), refine])
x_rates = [.035, .04, .045, .05, .055, .06, .065]
z_rates = [.025, .0275, .03, .0325, .035, .0375, .04]
sizes = [3, 5, 7, 9, 11]
shards = {}
for name, channel, ls, rates, seed in [
    ('production-x-small', 'x_only', [3, 5, 7], x_rates, 2026101571),
    ('production-x-L9', 'x_only', [9], x_rates, 2026101572),
    ('production-x-L11-a', 'x_only', [11], x_rates[:3], 2026101573),
    ('production-x-L11-b', 'x_only', [11], x_rates[3:5], 2026101574),
    ('production-x-L11-c', 'x_only', [11], x_rates[5:], 2026101575),
    ('production-z', 'z_only', sizes, z_rates, 2026101576),
]:
    shards[name] = dict(channel=channel, sizes=ls, rates=rates, shots=1000,
                        seed=seed, batch_size=50)
numerical = ('results/toric3d-h6-z-one-step-2026-10-15/analyze.py',
             'examples/fit_toric3d_h6_threshold.py',
             'examples/fit_toric3d_diagonal_feedback.py')
selection = dict(
    stage='production', pilot_excluded=True, frozen_before_sampling=True,
    frozen_at_unix=time.time(), bootstrap_replicates=200,
    channels={
        'x_only': dict(sizes=sizes, rates=x_rates, shots=1000, fit_rates=x_rates,
                       bootstrap_seed=2026101581),
        'z_only': dict(sizes=sizes, rates=z_rates, shots=1000, fit_rates=z_rates,
                       bootstrap_seed=2026101582)},
    shards=shards, source_hashes=sources,
    protocols={channel: m['protocol'] for channel, m in metas.items()},
    numerical_source_sha256={path: sha(ROOT / path) for path in numerical},
    pilot_input_sha256={f'pilot-{channel}': input_hashes(f'pilot-{channel}') for channel in 'xz'},
    additional_pilot_input_sha256={'pilot-z-refine': input_hashes('pilot-z-refine')},
    total_fault_records=70000, total_protocol_outcomes=140000,
    excluded_pilot_fault_records=2220,
    rationale='Freeze five matched sizes and 1000 independent records per point. '
              'Broad X pilot supports .035-.065 and preserves the preceding X grid. '
              'The 100-shot Z refinement places the crossing near .03-.035; .025-.04 '
              'centers seven rates around it, brackets size-order reversal and avoids '
              'the nearly saturated .045-.05 range. Six of seven Z rates match the '
              'preceding two-slice study. All 2220 pilot records are excluded from '
              'production estimates and fitting.',
    analysis_note='All numerical fit routines and guard policy are frozen before production; '
                  'input-validation and presentation code may finish independently. Existing '
                  'inverse-size corrections, drop-smallest and narrowed-window checks remain mandatory.')
original_path = STUDY / 'pre-optimization/selection.json'
if original_path.exists():
    original = json.loads(original_path.read_text())
    design_keys = ('channels', 'shards', 'bootstrap_replicates',
                   'total_fault_records', 'total_protocol_outcomes')
    assert all(selection[key] == original[key] for key in design_keys), \
        'A performance implementation change must not change the frozen statistical design.'
    assert selection['numerical_source_sha256'] == original['numerical_source_sha256']
    assert selection['protocols'] == original['protocols']
    selection['statistical_design_frozen_at_unix'] = original['frozen_at_unix']
    selection['original_selection_sha256'] = sha(original_path)
    selection['performance_refreeze_note'] = (
        'The original selection and all earlier outputs are archived under pre-optimization/. '
        'Only the implementation of exact binary matrix parity was optimized. The grid, '
        'seeds, sample counts, decoder settings, fit windows, bootstrap seeds and numerical '
        'policy are identical to the original frozen design. All pilots and production '
        'are rerun under the final source manifest; repeated streams are not pooled.')
(STUDY / 'selection.json').write_text(json.dumps(selection, indent=2) + '\n')
with tarfile.open(STUDY / 'verification/numerical-source.tar.gz', 'w:gz') as archive:
    for path in numerical:
        archive.add(ROOT / path, arcname=path, recursive=False)
base = ('env JULIA_CONDAPKG_BACKEND=Null '
        'JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python '
        'JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 '
        'PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/h6-full-prefix-mpl '
        'julia --compiled-modules=existing --project=. --startup-file=no '
        'results/toric3d-h6-full-prefix-2026-10-15/run_study.jl')
commands = {}
for name, config in shards.items():
    commands[name] = (base + ' --channel ' + config['channel']
        + ' --output results/toric3d-h6-full-prefix-2026-10-15/' + name
        + ' --sizes ' + ','.join(map(str, config['sizes']))
        + ' --rates ' + ','.join(map(str, config['rates']))
        + ' --shots 1000 --seed ' + str(config['seed']) + ' --batch-size 50')
(STUDY / 'verification/production-commands.json').write_text(json.dumps(commands, indent=2) + '\n')
(STUDY / 'verification/commands.txt').write_text(
    '# Frozen production commands; append --resume only for interrupted checkpoints.\n'
    '# Run from a project restored from verification/simulation-source.tar.gz.\n\n'
    + '\n\n'.join(commands.values()) + '\n')
print('Frozen 70-point, 70,000-record production selection before production sampling.')
