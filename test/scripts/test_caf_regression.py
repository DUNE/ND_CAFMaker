import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import unittest


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


if __name__ == "__main__":
    unittest.main()