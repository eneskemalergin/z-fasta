#!/usr/bin/env python3
"""Record the files used by stats verification for conservative resume checks."""

import json
import sys
from pathlib import Path

from build_identity import build_identity, file_identity


def verification_identity(root, binary, tools):
    root = Path(root)
    inputs = sorted((root / "tests/data").glob("*.fasta"))
    inputs += sorted((root / "bench/shared/cache/messy_fixtures").glob("*.fasta"))
    return {
        "subject": build_identity(root, binary),
        "inputs": {str(path): file_identity(path) for path in inputs},
        "tools": {path: file_identity(path) for path in tools},
    }


def can_resume(identity):
    if not isinstance(identity, dict):
        return False
    subject = identity.get("subject") or {}
    return bool(
        subject.get("checkout_revision")
        and subject.get("checkout_dirty") is False
        and subject.get("executable")
        and identity.get("inputs")
        and all(value is not None for value in identity["inputs"].values())
        and identity.get("tools")
        and all(value is not None for value in identity["tools"].values())
    )


if __name__ == "__main__":
    print(json.dumps(verification_identity(sys.argv[1], sys.argv[2], sys.argv[3:]), separators=(",", ":")))
