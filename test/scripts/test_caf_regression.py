import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import unittest

import ci_summary


RUNNER = Path(__file__).resolve().parents[2] / "scripts" / "caf-regression"
FILENAME = "MiniProdN5p1_NDComplex_FHC.0000965.CAF.root"
FAKE_APPTAINER = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import sys

bindings = {}
for position, argument in enumerate(sys.argv):
    if argument == "--bind":
        parts = sys.argv[position + 1].split(":")
        if len(parts) >= 2:
            bindings[parts[1]] = Path(parts[0])
work = bindings["/work"]
body = sys.argv[-2]
if "cmake -S" in body:
    stage = "build"
elif "makeCAF --fcl" in body:
    stage = "run"
else:
    stage = "verify"
with open(os.environ["MOCK_LOG"], "a") as output:
    output.write(stage + "\n")
if os.environ.get("MOCK_FAIL") == stage:
    sys.exit(23)
if os.environ.get("MOCK_SIGNAL"):
    os.kill(os.getppid(), int(os.environ["MOCK_SIGNAL"]))
    sys.exit(0)
if stage == "build":
    (work / "build").mkdir(exist_ok=True)
    (work / "install/bin").mkdir(parents=True, exist_ok=True)
    (work / "install/bin/ndcaf_setup.sh").touch()
elif stage == "run":
    assert "/work/inputs" in bindings
    (work / "outputs" / os.environ["MOCK_FILENAME"]).write_bytes(b"ROOT")
else:
    (work / "outputs" / (os.environ["MOCK_FILENAME"] + ".json")).write_text(
        os.environ.get("MOCK_JSON", '{"value": 1}\n')
    )
'''
FAKE_WGET = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import sys

with open(os.environ["MOCK_LOG"], "a") as output:
    output.write("download\n")
for argument in sys.argv:
    if argument.startswith("--directory-prefix="):
        Path(argument.split("=", 1)[1], "downloaded").touch()
if os.environ.get("MOCK_FAIL") == "download":
    sys.exit(23)
'''


class RegressionRunnerTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="caf runner test ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repo = self.root / "repo with spaces"
        (self.repo / "scripts").mkdir(parents=True)
        shutil.copy2(RUNNER, self.repo / "scripts/caf-regression")
        (self.repo / "test/references").mkdir(parents=True)
        self.reference = self.repo / "test/references" / (FILENAME + ".json")
        self.reference.write_text('{"value": 1}\n')
        self.reference.chmod(0o640)
        self.image = self.root / "image"
        self.image.mkdir()
        self.setup = self.root / "setup.sh"
        self.setup.touch()
        (self.repo / "test/regression-environment.conf").write_text(
            f"CAF_IMAGE='{self.image}'\nPYTHON_IMAGE='{self.image}'\n"
            f"CAF_DEPS_SETUP='{self.setup}'\nPYTHON_VERSION=3.15.2\n"
            f"CAF_INPUT_URL=https://example.invalid/\nCAF_FILENAME={FILENAME}\n"
        )
        self.bin = self.root / "bin"
        self.bin.mkdir()
        for name, content in (("apptainer", FAKE_APPTAINER), ("wget", FAKE_WGET)):
            executable = self.bin / name
            executable.write_text(content)
            executable.chmod(0o755)
        self.cache = self.root / "cache with spaces"
        self.cache.mkdir()
        (self.cache / "keep").write_text("cached data")
        self.scratch = self.root / "scratch"
        self.scratch.mkdir()
        self.log = self.root / "commands.log"
        self.environment = dict(os.environ)
        self.environment.pop("GITHUB_OUTPUT", None)
        self.environment.update(
            PATH=f"{self.bin}:{os.environ['PATH']}",
            TMPDIR=str(self.scratch),
            MOCK_LOG=str(self.log),
            MOCK_FILENAME=FILENAME,
        )

    def invoke(self, *arguments):
        return subprocess.run(
            ["bash", str(self.repo / "scripts/caf-regression"), *arguments],
            cwd=self.root,
            env=self.environment,
            capture_output=True,
            text=True,
            timeout=20,
        )

    def local_run(self, *arguments):
        return self.invoke(*arguments, "--cache-dir", str(self.cache))

    def assert_clean(self):
        self.assertEqual(list(self.scratch.iterdir()), [])
        self.assertEqual(list(self.reference.parent.iterdir()), [self.reference])
        self.assertEqual((self.cache / "keep").read_text(), "cached data")

    def test_check_reuses_cache_and_preserves_repository_artifacts(self):
        for name in ("build", "install", "outputs"):
            (self.repo / name).mkdir()
            (self.repo / name / "keep").write_text("untouched")
        result = self.local_run("--check")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("WARNING: Reusing", result.stderr)
        self.assertEqual(self.log.read_text(), "build\nrun\nverify\n")
        for name in ("build", "install", "outputs"):
            self.assertEqual((self.repo / name / "keep").read_text(), "untouched")
        self.assert_clean()

    def test_mismatch_never_updates_reference(self):
        self.environment["MOCK_JSON"] = '{"value": 2}\n'
        result = self.local_run("--check")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('+{"value": 2}', result.stdout)
        self.assertEqual(self.reference.read_text(), '{"value": 1}\n')
        self.assert_clean()

    def test_ci_failure_banner_follows_diff_and_is_absent_on_success(self):
        url = "https://github.com/DUNE/ND_CAFMaker/actions/runs/123#summary"
        self.environment["CI_SUMMARY_URL"] = url
        result = self.local_run("--check")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("VIEW THE SUMMARY", result.stdout)
        self.environment["MOCK_JSON"] = '{"value": 2}\n'
        result = self.local_run("--check")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("\033[1;31mREGRESSION FAILED - VIEW THE SUMMARY\033[0m", result.stdout)
        self.assertIn(url, result.stdout)
        self.assertLess(result.stdout.index('+{"value": 2}'), result.stdout.index("VIEW THE SUMMARY"))
        self.assertEqual(self.reference.read_text(), '{"value": 1}\n')
        self.assert_clean()

    def test_update_replaces_atomically_and_preserves_permissions(self):
        old_inode = self.reference.stat().st_ino
        self.environment["MOCK_JSON"] = '{"value": 2}\n'
        result = self.local_run("--update-reference")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.reference.read_text(), '{"value": 2}\n')
        self.assertNotEqual(self.reference.stat().st_ino, old_inode)
        self.assertEqual(self.reference.stat().st_mode & 0o777, 0o640)
        self.assert_clean()

    def test_unchanged_update_does_not_replace(self):
        old_inode = self.reference.stat().st_ino
        result = self.local_run("--update-reference")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.reference.stat().st_ino, old_inode)
        self.assert_clean()

    def test_failed_stage_preserves_reference_and_cleans_scratch(self):
        for stage in ("build", "run", "verify"):
            with self.subTest(stage=stage):
                self.environment["MOCK_FAIL"] = stage
                result = self.local_run("--update-reference")
                self.assertEqual(result.returncode, 23, result.stderr)
                self.assertEqual(self.reference.read_text(), '{"value": 1}\n')
                self.assert_clean()

    def test_new_cache_download_and_failed_download_is_retained(self):
        for failed in (False, True):
            with self.subTest(failed=failed):
                fresh = self.root / f"fresh cache {failed}"
                if failed:
                    self.environment["MOCK_FAIL"] = "download"
                result = self.invoke("--cache-dir", str(fresh))
                self.assertEqual(result.returncode, 23 if failed else 0, result.stderr)
                self.assertTrue((fresh / "downloaded").exists())
                self.assert_clean()

    def test_interruption_cleans_scratch_and_preserves_reference(self):
        for interrupt, expected_status in ((signal.SIGINT, 130), (signal.SIGTERM, 143)):
            with self.subTest(interrupt=interrupt):
                self.environment["MOCK_SIGNAL"] = str(int(interrupt))
                result = self.local_run("--update-reference")
                self.assertEqual(result.returncode, expected_status, result.stderr)
                self.assertEqual(self.reference.read_text(), '{"value": 1}\n')
                self.assert_clean()

    def test_failed_rename_removes_sibling_temporary_file(self):
        failing_rename = self.bin / "mv"
        failing_rename.write_text("#!/bin/bash\nexit 23\n")
        failing_rename.chmod(0o755)
        self.environment["MOCK_JSON"] = '{"value": 2}\n'
        result = self.local_run("--update-reference")
        self.assertEqual(result.returncode, 23, result.stderr)
        self.assertEqual(self.reference.read_text(), '{"value": 1}\n')
        self.assert_clean()

    def test_ci_stages_preserve_artifacts_and_diagnostics(self):
        work = self.root / "ci work"
        for stage in ("build", "run", "verify"):
            result = self.invoke(
                "--stage", stage, "--work-dir", str(work),
                "--cache-dir", str(self.cache),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assert_clean()
        self.assertTrue((work / "install/bin/ndcaf_setup.sh").exists())
        self.assertTrue((work / "outputs" / FILENAME).exists())
        self.assertTrue((work / "outputs/reference.diff").exists())
        self.environment["MOCK_JSON"] = '{"value": 2}\n'
        result = self.invoke("--stage", "verify", "--work-dir", str(work))
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('+{"value": 2}', (work / "outputs/reference.diff").read_text())
        self.assertEqual(self.reference.read_text(), '{"value": 1}\n')

    def test_ci_comparison_output_distinguishes_mismatch_from_failure(self):
        work = self.root / "verification work"
        (work / "outputs").mkdir(parents=True)
        (work / "outputs" / FILENAME).write_bytes(b"ROOT")
        github_output = self.root / "github-output"
        self.environment["GITHUB_OUTPUT"] = str(github_output)
        for scenario, expected, status in (("match", "match", 0), ("mismatch", "mismatch", 1), ("error", "error", 23)):
            with self.subTest(scenario=scenario):
                github_output.write_text("")
                self.environment.pop("MOCK_FAIL", None)
                self.environment.pop("MOCK_JSON", None)
                if scenario == "mismatch":
                    self.environment["MOCK_JSON"] = '{"value": 2}\n'
                elif scenario == "error":
                    self.environment["MOCK_FAIL"] = "verify"
                result = self.invoke("--stage", "verify", "--work-dir", str(work))
                self.assertEqual(result.returncode, status, result.stderr)
                self.assertEqual(github_output.read_text().splitlines()[-1], f"comparison={expected}")
                self.assert_clean()

    def test_invalid_cli_and_missing_cache_fail_before_running(self):
        invalid = (
            ("--stage", "verify", "--work-dir", str(self.root), "--update-reference"),
            ("--stage", "invalid", "--work-dir", str(self.root)),
            ("--cache-dir",),
            ("--unknown",),
        )
        for arguments in invalid:
            with self.subTest(arguments=arguments):
                self.assertNotEqual(self.invoke(*arguments).returncode, 0)
        self.assertFalse(self.log.exists())
        self.assert_clean()


class CISummaryTest(unittest.TestCase):
    def successful_needs(self, comparison="match"):
        needs = {stage: {"result": "success", "outputs": {}} for stage in ("resolve-inputs", "build", "run-cafmaker", "verify")}
        needs["prepare-pr-inputs"] = {"result": "skipped", "outputs": {}}
        needs["verify"]["outputs"]["comparison"] = comparison
        return needs

    def test_passed_summary(self):
        report = ci_summary.render_summary(self.successful_needs(), [])
        self.assertIn("## Passed", report)
        self.assertNotIn("--update-reference", report)

    def test_mismatch_summary_and_update_command(self):
        needs = self.successful_needs("mismatch")
        needs["verify"]["result"] = "failure"
        report = ci_summary.render_summary(needs, [], "REGRESSION TABLE")
        self.assertIn("## Regression Mismatch", report)
        self.assertIn("REGRESSION TABLE", report)
        self.assertIn("./scripts/caf-regression --update-reference", report)

    def test_upstream_failure_and_skips_are_not_regression_mismatches(self):
        needs = self.successful_needs("")
        needs["resolve-inputs"]["result"] = "failure"
        for stage in ("build", "run-cafmaker", "verify"):
            needs[stage]["result"] = "skipped"
        report = ci_summary.render_summary(needs, [])
        self.assertIn("## CI Failure", report)
        self.assertIn("Regression was not completed", report)
        self.assertNotIn("--update-reference", report)

    def test_dump_error_has_no_update_advice(self):
        needs = self.successful_needs("error")
        needs["verify"]["result"] = "failure"
        report = ci_summary.render_summary(needs, [])
        self.assertIn("Verification failed before", report)
        self.assertNotIn("--update-reference", report)

    def test_artifact_failure_alongside_mismatch_is_reported(self):
        needs = self.successful_needs("mismatch")
        needs["verify"]["result"] = "failure"
        jobs = [{"name": "Verify CAF output", "html_url": "https://github.com/example/repo/actions/runs/1/job/2", "steps": [
            {"name": "Dump and compare CAF output", "conclusion": "failure"},
            {"name": "Upload comparison diagnostics", "conclusion": "failure"},
        ]}]
        report = ci_summary.render_summary(needs, jobs)
        self.assertIn("## CI Failure", report)
        self.assertIn("Upload comparison diagnostics", report)
        self.assertIn("Intentional Reference Changes", report)

    def test_preparation_only_does_not_claim_regression_pass(self):
        needs = {stage: {"result": "skipped", "outputs": {}} for stage in ci_summary.STAGES}
        needs["prepare-pr-inputs"]["result"] = "success"
        report = ci_summary.render_summary(needs, [], preparation_only=True)
        self.assertIn("## Input Preparation Passed", report)
        self.assertNotIn("## Passed", report)
        self.assertIn("did not build or run", report)

    def test_cancelled_and_unknown_results_do_not_pass(self):
        for status in ("cancelled", "skipped", "unknown"):
            with self.subTest(status=status):
                needs = self.successful_needs()
                needs["build"]["result"] = status
                self.assertIn("## CI Failure", ci_summary.render_summary(needs, []))

    def test_structured_diff_add_remove_type_and_array_changes(self):
        report = ci_summary.render_regression(
            {"energy": 1, "removed": None, "vector": [1, 2], "type": True},
            {"energy": 2, "added": None, "vector": [1, 3, 4], "type": 1}, "diff",
        )
        self.assertIn("6 changed fields", report)
        self.assertIn("/vector/2", report)
        self.assertIn("Added", report)
        self.assertIn("Removed", report)
        self.assertIn("(missing)", report)

    def test_diff_html_is_escaped_and_bounded(self):
        reference = {str(position): 1 for position in range(60)}
        candidate = {str(position): "<script>|danger</script>" for position in range(60)}
        report = ci_summary.render_regression(reference, candidate, "<script>\n" * 250)
        self.assertIn("Showing the first 50 fields", report)
        self.assertIn("Preview truncated", report)
        self.assertNotIn("<script>", report)
        self.assertIn("&#124;", report)

    def test_formatting_only_mismatch(self):
        report = ci_summary.render_regression({"value": 1}, {"value": 1}, "diff")
        self.assertIn("formatting differences", report)

    def test_missing_comparison_and_reporting_failure(self):
        self.assertIn("## CI Failure", ci_summary.render_summary(self.successful_needs(""), []))
        report = ci_summary.render_summary(self.successful_needs(), [], reporting_errors=["Artifact download failed"])
        self.assertIn("## CI Failure", report)
        self.assertIn("Artifact download failed", report)


if __name__ == "__main__":
    unittest.main()