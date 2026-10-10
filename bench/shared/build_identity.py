#!/usr/bin/env python3
"""Record cheap checkout and executable facts; build facts require a declaration."""

import json
import os
import subprocess
import stat
import sys
from pathlib import Path


def git_value(root, *args):
    try:
        return subprocess.run(
            ["git", "-C", str(root), *args], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def file_identity(filename):
    try:
        path = Path(filename).resolve(strict=True)
        facts = path.stat()
        if not stat.S_ISREG(facts.st_mode):
            return None
        return {
            "path": str(path),
            "size": facts.st_size,
            "mtime_ns": facts.st_mtime_ns,
            "ctime_ns": facts.st_ctime_ns,
            "device": facts.st_dev,
            "inode": facts.st_ino,
        }
    except OSError:
        return None


def build_identity(root, binary):
    revision = git_value(root, "rev-parse", "HEAD")
    changes = git_value(root, "status", "--porcelain", "--untracked-files=normal")
    executable = file_identity(binary)
    return {
        "checkout_revision": revision,
        "checkout_dirty": None if changes is None else bool(changes),
        "executable": executable,
        "declared_build": {
            "revision": os.environ.get("BENCH_BUILD_REVISION") or None,
            "zig_version": os.environ.get("BENCH_BUILD_ZIG_VERSION") or None,
            "optimize": os.environ.get("BENCH_BUILD_OPTIMIZE") or None,
            "target": os.environ.get("BENCH_BUILD_TARGET") or None,
        },
    }


if __name__ == "__main__":
    print(json.dumps(build_identity(sys.argv[1], sys.argv[2]), separators=(",", ":")))
