import tempfile
import unittest
from pathlib import Path

from pipeline.phase2_build import provenance


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="phase2-provenance-", dir="/tmp/opencode")
        self.root = Path(self.temporary.name)
        self.work = self.root / "data/interim/scoring"
        self.work.mkdir(parents=True)

    def tearDown(self):
        self.temporary.cleanup()

    def test_manifest_archive_preserves_exact_whitespace_and_key_order(self):
        source = self.root / "v1.json"
        data = b'{ "z": 1, "a": [2,3] }\r\n'
        source.write_bytes(data)
        target = self.work / "band_manifest_v1.json"
        provenance.archive_manifest(source, target)
        self.assertEqual(target.read_bytes(), data)
        provenance.archive_manifest(source, target)
        self.assertEqual(target.read_bytes(), data)

    def test_differing_archive_is_not_overwritten(self):
        source, target = self.root / "v1.json", self.work / "archive.json"
        source.write_text("{}")
        target.write_text('{"user": true}')
        with self.assertRaisesRegex(ValueError, "differs"):
            provenance.archive_manifest(source, target)
        self.assertEqual(target.read_text(), '{"user": true}')

    def test_archive_in_place_is_rejected(self):
        source = self.root / "v1.json"
        source.write_text("{}")
        with self.assertRaisesRegex(ValueError, "in place"):
            provenance.archive_manifest(source, source)

    def test_only_descriptive_interim_outputs_are_allowed(self):
        for target in (self.root / "data/inputs/scoring", self.root / "data/processed", self.root / "data/interim", self.root):
            with self.subTest(path=target), self.assertRaisesRegex(ValueError, "descriptive"):
                provenance.check_workspace(target, self.root)

    def test_existing_unfamiliar_candidate_is_preserved(self):
        target = self.work / "common_stack_v2.tif"
        target.write_text("user work")
        with self.assertRaisesRegex(ValueError, "no known owner"):
            provenance.check_workspace(self.work, self.root)
        self.assertEqual(target.read_text(), "user work")

    def test_own_complete_and_partial_outputs_can_be_rerun(self):
        work, managed = provenance.check_workspace(self.work, self.root)
        provenance.begin_workspace(work, managed)
        (work / "common_stack_v2_source.tif").write_bytes(b"partial owned build")
        self.assertEqual(provenance.check_workspace(work, self.root), (work, managed))
        provenance.finish_workspace(work, managed)
        self.assertEqual(provenance.check_workspace(work, self.root), (work, managed))

    def test_changed_completed_output_is_not_discarded(self):
        work, managed = provenance.check_workspace(self.work, self.root)
        provenance.begin_workspace(work, managed)
        target = work / "score_build_report.md"
        target.write_text("generated")
        provenance.finish_workspace(work, managed)
        target.write_text("user review notes")
        with self.assertRaisesRegex(ValueError, "changed outside"):
            provenance.check_workspace(work, self.root)
        self.assertEqual(target.read_text(), "user review notes")

    def test_unrelated_files_are_not_claimed_or_deleted(self):
        unrelated = self.work / "approved_params.yaml"
        unrelated.write_text("user parameters")
        work, managed = provenance.check_workspace(self.work, self.root)
        provenance.begin_workspace(work, managed)
        provenance.finish_workspace(work, managed)
        self.assertNotIn(unrelated.name, managed)
        self.assertEqual(unrelated.read_text(), "user parameters")

    def test_symlink_output_and_symlink_directory_escape_are_rejected(self):
        outside = self.root / "outside"
        outside.mkdir()
        (self.work / "common_stack_v2.tif").symlink_to(outside / "victim")
        with self.assertRaisesRegex(ValueError, "symlink"):
            provenance.check_workspace(self.work, self.root)
        (self.work / "common_stack_v2.tif").unlink()
        (self.work / "reproducibility").symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "escaped"):
            provenance.check_workspace(self.work, self.root)

    def test_snapshot_files_are_checksum_identical(self):
        script = self.root / "pipeline/phase2_build/test.py"
        script.parent.mkdir(parents=True)
        script.write_text("# Synthetic test script\n")
        work, managed = provenance.check_workspace(self.work, self.root, code_files=(script,))
        provenance.begin_workspace(work, managed)
        provenance.snapshot_code((script,), work, self.root)
        target = work / "script_snapshot/pipeline/phase2_build/test.py"
        self.assertEqual(script.read_bytes(), target.read_bytes())

    def test_json_nan_duplicates_and_nonfinite_writes_are_rejected(self):
        path = self.root / "bad.json"
        for text in ('{"key":NaN}', '{"key":Infinity}', '{"key":1e400}', '{"key":1,"key":2}'):
            path.write_text(text)
            with self.assertRaises(ValueError):
                provenance.load_json(path)
        with self.assertRaises(ValueError):
            provenance.json_bytes({"key": float("nan")})

    def test_symlink_parent_within_workspace_is_not_followed(self):
        target = self.work / "unrelated-user-directory"
        target.mkdir()
        (self.work / "reproducibility").symlink_to(target, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink"):
            provenance.check_workspace(self.work, self.root)
        self.assertFalse(any(target.iterdir()))

    def test_fingerprints_are_stable_and_use_clear_relative_paths(self):
        path = self.root / "path with spaces.txt"
        path.write_text("synthetic data")
        observed = provenance.fingerprints([path, path], self.root)
        self.assertEqual(list(observed), [path.name])
        self.assertEqual(observed[path.name]["size_bytes"], 14)
        self.assertEqual(observed, provenance.fingerprints([path], self.root))
