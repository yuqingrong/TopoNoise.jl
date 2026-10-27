#!/usr/bin/env python3
"""Read-only exact-count audit of the archived GF(2) performance transition.

Run after all required new generations exist. JSON is written only to stdout;
redirect it to a verification report if desired. No streams are decoded,
sampled, pooled, resumed or modified by this program.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys
import tarfile
import tomllib


ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
ALLOWED_SOURCE_CHANGE = 'src/toric3d/decoder.jl'
BASELINE_ARCHIVE_FILES = 1364
SOURCE_FILES = 56
DESIGN_KEYS = (
    'stage', 'pilot_excluded', 'frozen_before_sampling', 'channels', 'shards',
    'bootstrap_replicates', 'protocols', 'numerical_source_sha256',
    'total_fault_records', 'total_protocol_outcomes', 'excluded_pilot_fault_records',
)
STATE_KEYS = ('completed_shots', 'logical_1_counts', 'logical_2_counts', 'failures', 'joint_counts')


class TransitionError(ValueError):
    pass


class Audit:
    def __init__(self):
        self.counts = Counter()
        self.inputs = {}

    def require(self, category, condition, message):
        if not condition:
            raise TransitionError(message)
        self.counts[category] += 1

    def read(self, path):
        path = Path(path).resolve()
        blob = path.read_bytes()
        value = digest(blob)
        self.require('read_only', path not in self.inputs or self.inputs[path] == value,
                     f'Input changed during audit: {path}')
        self.inputs[path] = value
        return blob

    def unchanged(self):
        for path, expected in self.inputs.items():
            self.require('read_only', digest(path.read_bytes()) == expected,
                         f'Input changed during audit: {path}')


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def load_json(audit, path):
    return json.loads(audit.read(path))


def compare_design(audit, original, current):
    for key in DESIGN_KEYS:
        audit.require('design', key in original and key in current and current[key] == original[key],
                      f'Frozen scientific design differs: {key}')
    for key in ('pilot_input_sha256', 'additional_pilot_input_sha256'):
        audit.require('design', isinstance(original.get(key), dict) and isinstance(current.get(key), dict)
                      and set(current[key]) == set(original[key]), f'Pilot dataset set differs: {key}')
    audit.require('design', current.get('statistical_design_frozen_at_unix') == original['frozen_at_unix'],
                  'New selection does not preserve the original statistical freeze time')


def compare_sources(audit, original, current):
    audit.require('sources', isinstance(original, dict) and isinstance(current, dict) and
                  len(original) == len(current) == SOURCE_FILES and set(original) == set(current),
                  'The source manifest must retain exactly the same 56 paths')
    changed = sorted(path for path in original if current[path] != original[path])
    audit.require('sources', changed == [ALLOWED_SOURCE_CHANGE],
                  f'Expected only {ALLOWED_SOURCE_CHANGE} to change; observed {changed}')
    return changed


def verify_source_archive(audit, path, manifest):
    # Read members without extracting files, and reject duplicate or extra paths.
    import io
    blob = audit.read(path)
    with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        audit.require('source_archive', all(member.isfile() for member in members) and
                      len(names) == len(set(names)) and set(names) == set(manifest),
                      f'Source archive member set differs: {path}')
        for member in members:
            audit.require('source_archive', digest(archive.extractfile(member).read()) == manifest[member.name],
                          f'Source archive bytes differ: {member.name}')
    return digest(blob)


def verify_baseline_archive(audit, baseline):
    manifest_path = baseline/'ARCHIVE.sha256.json'
    manifest = load_json(audit, manifest_path)
    files = {str(path.relative_to(baseline)) for path in baseline.rglob('*')
             if path.is_file() and path != manifest_path}
    audit.require('baseline_archive', isinstance(manifest, dict) and len(manifest) == BASELINE_ARCHIVE_FILES
                  and set(manifest) == files, 'Archived baseline is not the complete 1,364-file manifest')
    for relative, expected in manifest.items():
        audit.require('baseline_archive', digest(audit.read(baseline/relative)) == expected,
                      f'Archived baseline file changed: {relative}')
    return digest(audit.read(manifest_path))


def read_checkpoint(audit, path):
    blob = audit.read(path)
    checkpoint = tomllib.loads(blob.decode())
    audit.require('checkpoint', checkpoint.get('format_version') == 1 and
                  type(checkpoint.get('generation')) is int and checkpoint['generation'] >= 0,
                  f'Invalid checkpoint format or generation: {path}')
    payload, occurrences = re.subn(rb'^payload_sha256 = "[0-9a-f]{64}"\n', b'', blob, flags=re.M)
    audit.require('checkpoint', occurrences == 1 and digest(payload) == checkpoint.get('payload_sha256'),
                  f'Checkpoint payload checksum differs: {path}')
    config = checkpoint['config']
    size_count = len(config['sizes'])*len(config['rates'])
    audit.require('checkpoint', len(checkpoint['states']) == len(checkpoint['point_seeds_hex']) == size_count,
                  f'Checkpoint grid/state/seed dimensions differ: {path}')
    audit.require('checkpoint', len(set(checkpoint['point_seeds_hex'])) == size_count and
                  all(re.fullmatch('[0-9a-f]{16}', seed) for seed in checkpoint['point_seeds_hex']),
                  f'Invalid or repeated exact point seed: {path}')
    for index, state in enumerate(checkpoint['states']):
        shots = state['completed_shots']
        audit.require('checkpoint', type(shots) is int and 0 <= shots <= config['shots'],
                      f'Invalid completed record count: {path}, point {index+1}')
        for key, length in (('logical_1_counts', 2), ('logical_2_counts', 2), ('failures', 2), ('joint_counts', 4)):
            values = state[key]
            audit.require('checkpoint', len(values) == length and all(type(value) is int and 0 <= value <= shots for value in values),
                          f'Invalid {key}: {path}, point {index+1}')
        audit.require('checkpoint', sum(state['joint_counts']) == shots, f'Joint total differs: {path}, point {index+1}')
        for arm in range(2):
            one, two, failure = state['logical_1_counts'][arm], state['logical_2_counts'][arm], state['failures'][arm]
            audit.require('checkpoint', max(one,two) <= failure <= min(shots,one+two) and
                          failure == sum(count for mask,count in enumerate(state['joint_counts']) if mask&(1<<arm)),
                          f'Logical/joint marginal differs: {path}, point {index+1}, arm {arm}')
    return checkpoint, blob


def snapshot_name(checkpoint):
    return f"checkpoint-{checkpoint['generation']:06d}-{checkpoint['payload_sha256'][:12]}.toml"


def find_matching_generation(audit, directory, generation):
    matches = sorted((directory/'checkpoints').glob(f'checkpoint-{generation:06d}-*.toml'))
    audit.require('generation', len(matches) == 1,
                  f'{directory.name}: expected exactly one new snapshot at generation {generation}; '
                  f'found {len(matches)}. The required generation must be completed before this audit passes.')
    return matches[0]


def compare_scientific_checkpoint(audit, name, old, new):
    for key in ('generation', 'config', 'point_seeds_hex', 'protocol'):
        audit.require('scientific_counts', new[key] == old[key], f'{name}: checkpoint {key} differs')
    audit.require('scientific_counts', len(new['states']) == len(old['states']), f'{name}: state count differs')
    for index, (before, after) in enumerate(zip(old['states'],new['states'])):
        for key in STATE_KEYS:
            audit.require('scientific_counts', after[key] == before[key],
                          f'{name}: point {index+1} differs in {key} at generation {old["generation"]}')


def check_transition(study=STUDY, baseline=None):
    study = Path(study).resolve()
    baseline = Path(baseline).resolve() if baseline else study/'pre-optimization'
    audit = Audit()
    baseline_archive_sha256 = verify_baseline_archive(audit, baseline)
    original = load_json(audit, baseline/'selection.json')
    current = load_json(audit, study/'selection.json')
    compare_design(audit, original, current)
    original_sha256 = digest(audit.read(baseline/'selection.json'))
    audit.require('design', current.get('original_selection_sha256') == original_sha256,
                  'New selection does not bind the exact archived original selection')
    changed = compare_sources(audit, original['source_hashes'], current['source_hashes'])
    source_archives = {}
    for label, directory, selection in (('baseline',baseline,original),('current',study,current)):
        manifest = load_json(audit, directory/'verification/simulation-source-sha256.json')
        audit.require('sources', manifest == selection['source_hashes'], f'{label}: selection/source manifest differs')
        source_archives[label] = verify_source_archive(audit, directory/'verification/simulation-source.tar.gz', manifest)
    for path, expected in current['source_hashes'].items():
        audit.require('sources', digest(audit.read(ROOT/path)) == expected, f'Current source differs from final manifest: {path}')
    for path, expected in current['numerical_source_sha256'].items():
        audit.require('sources', digest(audit.read(ROOT/path)) == expected, f'Frozen numerical source changed: {path}')

    pilot_names = set(original['pilot_input_sha256']) | set(original['additional_pilot_input_sha256'])
    production_names = set(original['shards'])
    expected_names = pilot_names | production_names
    authoritative = {path.parent.name:path for path in baseline.glob('*/checkpoint.toml')}
    audit.require('checkpoint_set', set(authoritative) == expected_names and len(expected_names) == 9,
                  'Must compare every archived pilot and production authoritative checkpoint (nine datasets)')
    comparisons = []
    for name in sorted(expected_names):
        before, before_bytes = read_checkpoint(audit, authoritative[name])
        baseline_snapshot = baseline/name/'checkpoints'/snapshot_name(before)
        audit.require('checkpoint', audit.read(baseline_snapshot) == before_bytes,
                      f'{name}: archived authoritative checkpoint differs from its immutable snapshot')
        new_path = find_matching_generation(audit, study/name, before['generation'])
        after, after_bytes = read_checkpoint(audit, new_path)
        audit.require('checkpoint', new_path.name == snapshot_name(after), f'{name}: new snapshot filename/checksum differs')
        audit.require('sources', before['source_hashes'] == original['source_hashes'] and
                      after['source_hashes'] == current['source_hashes'], f'{name}: snapshot uses the wrong transition source')
        compare_scientific_checkpoint(audit, name, before, after)
        comparisons.append(dict(dataset=name, category='pilot' if name in pilot_names else 'production',
                                channel=before['config']['channel'], generation=before['generation'],
                                grid_points=len(before['states']), completed_fault_records=sum(s['completed_shots'] for s in before['states']),
                                baseline_checkpoint=str(authoritative[name]), baseline_checkpoint_sha256=digest(before_bytes),
                                final_source_snapshot=str(new_path), final_source_snapshot_sha256=digest(after_bytes),
                                compared_state_fields=list(STATE_KEYS), exact_match=True))
    audit.unchanged()
    return dict(status='passed', verification='same-generation exact scientific-count reproduction after GF(2) parity optimization',
                baseline_directory=str(baseline), current_directory=str(study),
                baseline_archive_sha256=baseline_archive_sha256, original_selection_sha256=original_sha256,
                final_selection_sha256=digest(audit.read(study/'selection.json')),
                scientific_design_keys=list(DESIGN_KEYS), changed_source_paths=changed,
                source_archive_sha256=source_archives, comparisons=comparisons,
                reproduced_pilot_records=sum(x['completed_fault_records'] for x in comparisons if x['category']=='pilot'),
                reproduced_production_records=sum(x['completed_fault_records'] for x in comparisons if x['category']=='production'),
                assertions=sum(audit.counts.values()), assertions_by_category=dict(audit.counts),
                ignored_for_scientific_count_comparison=['elapsed_seconds','runtime','timestamps','source_hashes','decoder_build_timing'],
                source_validation='Source hashes are bound separately to original/final manifests; exactly decoder.jl may change.',
                independence_limit='Compares saved scientific counts at equal generations; does not establish equality of every individual recovery.',
                pooling='none; archived baseline counts are verification evidence only')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    try:
        report = check_transition(args.study,args.baseline)
    except (TransitionError, OSError, KeyError, ValueError, tarfile.TarError) as error:
        print(json.dumps(dict(status='failed', error=str(error))), file=sys.stderr)
        return 1
    print(json.dumps(report,indent=2,sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
