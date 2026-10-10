#!/usr/bin/env bash
# Consumer checks use actual process status and parsed reports together.
set -euo pipefail
binary=${1:-zig-out/bin/z-fasta}
python3 - "$binary" <<'PY'
import json
import subprocess
import sys
import tempfile
from pathlib import Path

binary = sys.argv[1]
for command in ("index", "get", "stats", "validate"):
    result = subprocess.run([binary, command, "--error-format=json"], capture_output=True, check=False)
    assert result.returncode == 1 and result.stdout == b""
    report = json.loads(result.stderr)
    assert report["schema_version"] == "v2"
    assert report["operation"]["name"] == command
    assert report["operation"]["complete"] is False
    assert report["accepted"] is False
    assert report["failure"]["primary"]["code"] == "usage.missing_argument"
    # Missing/truncated machine objects cannot be accepted just from status.
    for broken in (b"", result.stderr[:-2]):
        try:
            json.loads(broken)
        except (ValueError, UnicodeError):
            pass
        else:
            raise AssertionError("Incomplete JSON was accepted")

result = subprocess.run([binary, "get", "tests/data/simple.fasta", "seq1:0-1", "--error-format=json"], capture_output=True, check=False)
assert result.returncode == 1 and result.stdout == b""
report = json.loads(result.stderr)
assert report["failure"]["primary"]["code"] == "request.start_before_one"
assert report["failure"]["primary"]["data"]["start"] == "0"

# Warning acceptance is explicit. Validator data remains a separate v1 report.
with tempfile.TemporaryDirectory(prefix="zfasta-diagnostic-consumer-") as directory:
    warning = Path(directory) / "warning.fa"
    warning.write_bytes(b">empty\n")
    result = subprocess.run([binary, "validate", str(warning), "--json", "--summary", "--error-format=json"], capture_output=True, check=False)
    assert result.returncode == 2 and result.stderr == b""
    assert json.loads(result.stdout)["schema_version"] == "v1"
print("diagnostic consumer: passed")
PY
