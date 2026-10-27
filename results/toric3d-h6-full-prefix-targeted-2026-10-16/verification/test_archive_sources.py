"""Archive-gate regressions with temporary fixture files; no real data writes."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

STUDY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('archive_sources_under_test', STUDY / 'archive_sources.py')
ARCHIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ARCHIVE)


class ArchiveGateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='h6-archive-gate-test-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.study = self.root / 'results/study'
        self.verification = self.study / 'verification'
        self.verification.mkdir(parents=True)
        self.analysis = self.study / 'analysis'
        self.analysis.mkdir()
        self.shard = self.study / 'production-fixture'
        self.shard.mkdir()
        for name in ('raw.csv', 'checkpoint.toml', 'metadata.toml'):
            (self.shard / name).write_text('retained fixture\n')
        for name in ('check_production.py', 'check_numerics.py'):
            (self.verification / name).write_text('# retained checker fixture\n')
        self.source = self.root / 'source.py'
        self.source.write_text('# retained numerical fixture\n')
        (self.study / 'selection.json').write_text(json.dumps({'jobs': [{'name': 'production-fixture'}]}))
        (self.analysis / 'input.json').write_text('{"fixture": true}\n')
        self.report = {'status': 'complete_with_numerical_failures', 'completed_at_unix': 1,
                       'provenance': {'numerical_sources_sha256': {'source.py': ARCHIVE.sha(self.source)},
                                      'input_sha256': ARCHIVE.sha(self.analysis / 'input.json')}}
        (self.analysis / 'report.json').write_text(json.dumps(self.report))
        common = {'status': 'passed', 'selection_sha256': ARCHIVE.sha(self.study / 'selection.json')}
        self.production = {**common, 'checker_sha256': ARCHIVE.sha(self.verification / 'check_production.py'),
                           'shards': [{'name': 'production-fixture',
                                       'raw_sha256': ARCHIVE.sha(self.shard / 'raw.csv'),
                                       'checkpoint_sha256': ARCHIVE.sha(self.shard / 'checkpoint.toml'),
                                       'metadata_sha256': ARCHIVE.sha(self.shard / 'metadata.toml')}]}
        self.numerical = {**common, 'report_sha256': ARCHIVE.sha(self.analysis / 'report.json'),
                          'input_sha256': ARCHIVE.sha(self.analysis / 'input.json'),
                          'checker_sha256': ARCHIVE.sha(self.verification / 'check_numerics.py')}
        self.write_audits()
        self.root_patch = patch.object(ARCHIVE, 'ROOT', self.root)
        self.study_patch = patch.object(ARCHIVE, 'HERE', self.study)
        self.root_patch.start(); self.study_patch.start()
        self.addCleanup(self.root_patch.stop); self.addCleanup(self.study_patch.stop)

    def write_audits(self):
        (self.verification / 'production-audit.json').write_text(json.dumps(self.production))
        (self.verification / 'numerical-audit.json').write_text(json.dumps(self.numerical))

    def test_accounted_numerical_failures_allow_source_archiving(self):
        result = ARCHIVE.require_final_audits()
        self.assertEqual(result['final_report_sha256'], ARCHIVE.sha(self.analysis / 'report.json'))

    def test_each_production_artifact_change_invalidates_audit(self):
        for name in ('raw.csv', 'checkpoint.toml', 'metadata.toml'):
            with self.subTest(name=name):
                path = self.shard / name
                original = path.read_bytes()
                path.write_bytes(original + b'changed after audit\n')
                with self.assertRaisesRegex(ValueError, 'Production audit is stale'):
                    ARCHIVE.require_final_audits()
                path.write_bytes(original)

    def test_missing_or_duplicate_audited_shard_is_rejected(self):
        original = list(self.production['shards'])
        for shards in ([], original * 2):
            self.production['shards'] = shards; self.write_audits()
            with self.assertRaisesRegex(ValueError, 'every planned shard exactly once'):
                ARCHIVE.require_final_audits()

    def test_optional_production_checker_binding_is_enforced(self):
        path = self.verification / 'check_production.py'
        path.write_text('# changed checker\n')
        with self.assertRaisesRegex(ValueError, 'Production audit checker changed'):
            ARCHIVE.require_final_audits()
        self.production.pop('checker_sha256'); self.write_audits()
        ARCHIVE.require_final_audits()

    def test_numerical_input_change_is_rejected(self):
        (self.analysis / 'input.json').write_text('{"fixture": false}\n')
        with self.assertRaisesRegex(ValueError, 'Numerical input binding is stale'):
            ARCHIVE.require_final_audits()

    def test_numerical_source_and_checker_changes_are_rejected(self):
        self.source.write_text('# changed numerical source\n')
        with self.assertRaisesRegex(ValueError, 'Numerical source changed'):
            ARCHIVE.require_final_audits()
        self.source.write_text('# retained numerical fixture\n')
        (self.verification / 'check_numerics.py').write_text('# changed checker\n')
        with self.assertRaisesRegex(ValueError, 'Numerical audit checker changed'):
            ARCHIVE.require_final_audits()

    def test_running_report_is_rejected(self):
        self.report['status'] = 'running'
        (self.analysis / 'report.json').write_text(json.dumps(self.report))
        with self.assertRaisesRegex(ValueError, 'Analysis is not complete'):
            ARCHIVE.require_final_audits()

    def test_unsafe_paths_are_rejected(self):
        for value in ('../outside', '/outside', 'dir/../../outside', '.'):
            with self.subTest(path=value), self.assertRaises(ValueError):
                ARCHIVE.checked_relative(value)


if __name__ == '__main__':
    unittest.main(verbosity=2)
