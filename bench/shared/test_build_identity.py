#!/usr/bin/env python3
"""Check build metadata without building or measuring a command."""

import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from build_identity import build_identity, file_identity


class BuildIdentityTests(unittest.TestCase):
    def git(self, root, *args):
        return subprocess.run(
            ["git", "-c", "user.name=Metadata test", "-c", "user.email=test@example.invalid",
             "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", "-C", str(root), *args],
            check=True, capture_output=True, text=True,
        ).stdout.strip()

    def test_file_facts_reject_directories_and_keep_regular_symlinks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "subject"
            binary.write_bytes(b"fixture")
            link = root / "link"
            link.symlink_to(binary)
            broken = root / "broken"
            broken.symlink_to(root / "missing")
            for path in (root, "", root / "missing", broken):
                with self.subTest(path=str(path)):
                    self.assertIsNone(file_identity(path))
            self.assertEqual(file_identity(binary), file_identity(link))
            self.assertEqual(file_identity(binary)["size"], 7)

    def test_revisions_dirty_state_and_unknown_build_are_distinct(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "subject"
            binary.write_bytes(b"fixture")
            self.git(root, "init", "-q")
            self.git(root, "add", "subject")
            self.git(root, "commit", "-qm", "test: first revision")
            with patch.dict(os.environ, {key: "" for key in (
                "BENCH_BUILD_REVISION", "BENCH_BUILD_ZIG_VERSION", "BENCH_BUILD_OPTIMIZE", "BENCH_BUILD_TARGET"
            )}):
                first = build_identity(root, binary)
                self.assertFalse(first["checkout_dirty"])
                self.assertTrue(all(value is None for value in first["declared_build"].values()))
                binary.write_bytes(b"changed fixture")
                dirty = build_identity(root, binary)
                self.assertTrue(dirty["checkout_dirty"])
                self.assertNotEqual(first["executable"], dirty["executable"])
                self.git(root, "add", "subject")
                self.git(root, "commit", "-qm", "test: second revision")
                second = build_identity(root, binary)
                self.assertNotEqual(first["checkout_revision"], second["checkout_revision"])
                self.assertFalse(second["checkout_dirty"])
                with tempfile.TemporaryDirectory() as unknown:
                    missing = build_identity(unknown, Path(unknown) / "missing")
                    self.assertIsNone(missing["checkout_revision"])
                    self.assertIsNone(missing["checkout_dirty"])
                    self.assertIsNone(missing["executable"])


if __name__ == "__main__":
    unittest.main()
