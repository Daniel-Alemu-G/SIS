import argparse
import contextlib
import io
import unittest
from unittest import mock

from pipeline.phase2_build import build_scores, provenance, technical
from pipeline.phase2_build.tests import test_build_scores as build_fixtures
from pipeline.phase2_promote import finalize_scores


class FinalizationTests(unittest.TestCase):
    def setUp(self):
        # Reuse the tiny accepted-v1/config fixtures, not the real project's
        # data or dimensions. Composition avoids inheriting/rerunning its tests.
        self.fixture = build_fixtures.CandidateBuildTests("test_complete_candidate_build_preserves_v1_and_never_promotes")
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        self.root, self.candidate_dir = self.fixture.root, self.fixture.args.work_dir
        self.report = self.fixture.build()
        self.candidate = provenance.load_json(self.candidate_dir / "band_manifest.json")
        self.original_v1 = self.fixture.original_manifest
        self.patch = mock.patch.object(finalize_scores, "ROOT", self.root)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        for name in ("__init__.py", "finalize_scores.py"):
            path = self.root / "pipeline/phase2_promote" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# Synthetic publisher provenance fixture\n")
        self.work_dir = self.root / "data/interim/score_promote"
        self.work_dir.mkdir(parents=True)
        self.approval_path = self.work_dir / "owner_approval.json"
        self.approval = {
            "approved": True, "authority": "Synthetic project owner", "date": "2026-11-25",
            "statement": "Approve this exact advisory candidate with the disclosed flags; finalize Phase 2.",
            "scope": finalize_scores.SCOPE, "candidate_raster_sha256": self.candidate["raster_sha256"],
            "candidate_manifest_sha256": provenance.sha256(self.candidate_dir / "band_manifest.json"),
            "accepted_review_flags": [row["code"] for row in self.report["diagnostics"]["warnings"]],
        }
        provenance.write_json(self.approval_path, self.approval)
        self.args = argparse.Namespace(candidate_dir=self.candidate_dir, work_dir=self.work_dir, approval=self.approval_path)
        self.processed = self.root / "data/processed"
        self.candidate_before = self.files(self.candidate_dir)

    def files(self, directory=None):
        return {path.relative_to(self.root).as_posix(): (path.stat().st_size, path.stat().st_mtime_ns, provenance.sha256(path))
                for path in (directory or self.root).rglob("*") if path.is_file()}

    def finalize(self):
        with mock.patch.object(build_scores, "code_files", return_value=self.fixture.scripts), contextlib.redirect_stdout(io.StringIO()):
            return finalize_scores.finalize(self.args)

    def validate(self):
        with mock.patch.object(build_scores, "code_files", return_value=self.fixture.scripts):
            return finalize_scores.validate_only(self.processed)

    def test_complete_finalization_archives_v1_and_preserves_exact_approved_raster(self):
        result = self.finalize()
        self.assertEqual(result["status"], "finalized")
        self.assertEqual(result["band_count"], 52)
        manifest = provenance.load_json(self.processed / "band_manifest.json")
        self.assertEqual(manifest["status"], "accepted")
        self.assertEqual(manifest["bands"], self.candidate["bands"])
        self.assertEqual(manifest["scoring"], self.candidate["scoring"])
        self.assertEqual(manifest["build_inputs"], self.candidate["build_inputs"])
        self.assertEqual(manifest["raster_file"], "data/processed/common_stack_v2.tif")
        self.assertEqual(manifest["parent_stack"]["manifest_file"], "data/processed/band_manifest_v1.json")
        self.assertEqual((self.processed / "band_manifest_v1.json").read_bytes(), self.original_v1)
        self.assertEqual((self.processed / "band_manifest_v2.json").read_bytes(), (self.processed / "band_manifest.json").read_bytes())
        self.assertEqual(provenance.sha256(self.processed / "common_stack_v2.tif"), self.candidate["raster_sha256"])
        self.assertEqual(provenance.sha256(self.fixture.args.source_stack), self.fixture.v1["raster_sha256"])
        self.assertEqual(self.candidate_before, self.files(self.candidate_dir))
        self.assertEqual((self.processed / "stack_v2_approval.json").read_bytes(), self.approval_path.read_bytes())
        report = provenance.load_json(self.processed / "score_build_report.json")
        self.assertEqual(report["status"], "accepted_owner_review")
        self.assertEqual(report["diagnostics"], self.report["diagnostics"])
        self.assertEqual(report["reproducibility"], self.report["reproducibility"])
        self.assertIn("accepted_owner_review", (self.processed / "score_build_report.md").read_text())
        self.assertNotIn("candidate only; not promoted", (self.processed / "score_build_report.md").read_text())
        self.assertEqual(provenance.load_json(self.work_dir / finalize_scores.STATE_FILE)["status"], "accepted")

    def test_archive_exists_and_shared_is_still_v1_during_full_post_copy_acceptance(self):
        observed = []
        original = technical.validate_package

        def validate(*args, **kwargs):
            observed.append((self.processed / "band_manifest.json").read_bytes())
            self.assertEqual((self.processed / "band_manifest_v1.json").read_bytes(), self.original_v1)
            self.assertTrue((self.processed / "common_stack_v2.tif").is_file())
            return original(*args, **kwargs)

        with mock.patch.object(technical, "validate_package", side_effect=validate):
            self.finalize()
        self.assertEqual(observed, [self.original_v1])

    def test_published_package_revalidation_is_read_only(self):
        self.finalize()
        before = self.files()
        result = self.validate()
        self.assertEqual(result["status"], "passed_read_only")
        self.assertGreaterEqual(result["scalar_checks"], 32)
        self.assertEqual(result["review_flags"], self.approval["accepted_review_flags"])
        self.assertEqual(before, self.files())

    def test_completed_rerun_is_noop_and_byte_identical(self):
        self.finalize()
        before = self.files()
        self.assertEqual(self.finalize()["status"], "already_finalized")
        self.assertEqual(before, self.files())

    def test_approval_false_or_missing_identity_fails_before_outputs_are_claimed(self):
        for key, value in (("approved", False), ("candidate_raster_sha256", "0" * 64),
                           ("candidate_manifest_sha256", "1" * 64), ("authority", ""),
                           ("date", "invalid"), ("scope", "Start Phase 3")):
            with self.subTest(key=key):
                changed = {**self.approval, key: value}
                provenance.write_json(self.approval_path, changed)
                with self.assertRaises(ValueError):
                    self.finalize()
                self.assertFalse((self.work_dir / finalize_scores.STATE_FILE).exists())
                self.assertEqual((self.processed / "band_manifest.json").read_bytes(), self.original_v1)

    def test_every_scientific_flag_requires_explicit_acceptance(self):
        self.assertTrue(self.approval["accepted_review_flags"])
        changed = {**self.approval, "accepted_review_flags": self.approval["accepted_review_flags"][:-1]}
        provenance.write_json(self.approval_path, changed)
        with self.assertRaisesRegex(ValueError, "every disclosed review flag"):
            self.finalize()
        self.assertFalse((self.processed / "common_stack_v2.tif").exists())

    def test_changed_candidate_report_is_not_promoted(self):
        path = self.candidate_dir / "score_build_report.json"
        path.write_text(path.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "changed outside"):
            self.finalize()
        self.assertEqual((self.processed / "band_manifest.json").read_bytes(), self.original_v1)

    def test_incomplete_build_is_not_promoted(self):
        state_path = self.candidate_dir / provenance.STATE_FILE
        state = provenance.load_json(state_path)
        state["status"] = "building"
        provenance.write_json(state_path, state)
        with self.assertRaisesRegex(ValueError, "completed candidate"):
            self.finalize()

    def test_changed_threshold_input_is_not_silently_used_or_published(self):
        self.fixture.args.params.write_text(self.fixture.args.params.read_text() + "# New owner settings\n")
        with self.assertRaisesRegex(ValueError, "input preservation"):
            self.finalize()
        self.assertEqual((self.processed / "band_manifest.json").read_bytes(), self.original_v1)

    def test_unfamiliar_existing_outputs_are_not_overwritten(self):
        path = self.processed / "common_stack_v2.tif"
        path.write_bytes(b"Unrelated user raster")
        with self.assertRaisesRegex(ValueError, "no known owner"):
            self.finalize()
        self.assertEqual(path.read_bytes(), b"Unrelated user raster")
        self.assertFalse((self.work_dir / finalize_scores.STATE_FILE).exists())

    def test_matching_existing_v1_archive_is_preserved_without_rewriting(self):
        path = self.processed / "band_manifest_v1.json"
        path.write_bytes(self.original_v1)
        before = path.stat().st_mtime_ns
        self.finalize()
        self.assertEqual(path.stat().st_mtime_ns, before)

    def test_conflicting_v1_archive_is_never_overwritten(self):
        path = self.processed / "band_manifest_v1.json"
        path.write_bytes(b"Other historical manifest")
        with self.assertRaisesRegex(ValueError, "versioned v1 manifest differs"):
            self.finalize()
        self.assertEqual(path.read_bytes(), b"Other historical manifest")
        self.assertFalse((self.work_dir / finalize_scores.STATE_FILE).exists())

    def test_unfamiliar_shared_manifest_is_never_reverted(self):
        path = self.processed / "band_manifest.json"
        path.write_bytes(b"Newer user publication")
        with self.assertRaisesRegex(ValueError, "neither the immutable v1"):
            self.finalize()
        self.assertEqual(path.read_bytes(), b"Newer user publication")

    def test_symlink_and_atomic_scratch_outputs_are_rejected(self):
        victim = self.root / "private.txt"
        victim.write_bytes(b"Unrelated source")
        for name in ("common_stack_v2.tif", "band_manifest.json.tmp"):
            with self.subTest(name=name):
                target = self.processed / name
                target.symlink_to(victim)
                with self.assertRaisesRegex(ValueError, "symlink|outside the project"):
                    self.finalize()
                self.assertEqual(victim.read_bytes(), b"Unrelated source")
                target.unlink()

    def test_workspace_cannot_escape_interim_or_overlap_candidate(self):
        for path in (self.root / "data/processed/publish", self.candidate_dir / "publish",
                     self.root / "data/interim/../../data/inputs/publish"):
            with self.subTest(path=path):
                self.args.work_dir = path
                with self.assertRaisesRegex(ValueError, "separate descriptive directory"):
                    self.finalize()
        self.assertFalse((self.processed / "common_stack_v2.tif").exists())

    def test_failed_full_validation_does_not_advance_shared_manifest_and_can_resume(self):
        with mock.patch.object(technical, "validate_package", side_effect=RuntimeError("Synthetic interruption")):
            with self.assertRaisesRegex(RuntimeError, "Synthetic interruption"):
                self.finalize()
        self.assertEqual((self.processed / "band_manifest.json").read_bytes(), self.original_v1)
        self.assertEqual((self.processed / "band_manifest_v1.json").read_bytes(), self.original_v1)
        self.assertEqual(provenance.load_json(self.work_dir / finalize_scores.STATE_FILE)["status"], "finalizing")
        self.assertEqual(self.finalize()["status"], "finalized")

    def test_externally_changed_owned_output_is_not_overwritten_on_resume(self):
        with mock.patch.object(technical, "validate_package", side_effect=RuntimeError("Synthetic interruption")):
            with self.assertRaises(RuntimeError):
                self.finalize()
        path = self.processed / "common_stack_v2.tif"
        path.write_bytes(b"User replaced interrupted output")
        with self.assertRaisesRegex(ValueError, "changed outside"):
            self.finalize()
        self.assertEqual(path.read_bytes(), b"User replaced interrupted output")

    def test_corrupt_published_raster_fails_read_only_without_regeneration(self):
        self.finalize()
        path = self.processed / "common_stack_v2.tif"
        path.write_bytes(b"Corrupt externally modified raster")
        before = self.files()
        with self.assertRaisesRegex(ValueError, "published deliverables"):
            self.validate()
        self.assertEqual(before, self.files())

    def test_published_diagnostic_changes_are_rejected(self):
        self.finalize()
        path = self.processed / "score_build_report.json"
        changed = provenance.load_json(path)
        changed["diagnostics"]["warnings"] = []
        provenance.write_json(path, changed)
        with self.assertRaisesRegex(ValueError, "published deliverables"):
            self.validate()

    def test_candidate_archive_and_lookup_remain_sufficient_for_old_v1_contract(self):
        self.finalize()
        accepted = provenance.load_json(self.processed / "band_manifest.json")
        from config import stack, stack_v2
        archived = provenance.load_json(self.processed / "band_manifest_v1.json")
        stack.validate_manifest(archived)
        stack_v2.validate_manifest(accepted, v1_manifest=archived)
        self.assertEqual(accepted["bands"][:len(stack.CORE_BAND_NAMES)], archived["bands"])
        self.assertEqual(accepted["lookup"], archived["lookup"])


if __name__ == "__main__":
    unittest.main()
