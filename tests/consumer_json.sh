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
    assert "schema_version" not in report and report["catalog_revision"] == "0.4.0-1"
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

# Warning acceptance is explicit: the validate report on stdout says warnings, and the status is 2.
with tempfile.TemporaryDirectory(prefix="zfasta-diagnostic-consumer-") as directory:
    reference = Path(directory) / "reference.fa"
    reference.write_bytes(b">seq1\nACGT\n")
    Path(str(reference) + ".fai").write_bytes(b"seq1\t4\t6\t4\t5\n")
    result = subprocess.run([binary, "get", str(reference), "seq1:0-1", "--error-format=json"], capture_output=True, check=False)
    assert result.returncode == 1 and result.stdout == b""
    report = json.loads(result.stderr)
    assert report["failure"]["primary"]["code"] == "request.start_before_one"
    assert report["failure"]["primary"]["data"]["start"] == "0"

    warning = Path(directory) / "warning.fa"
    warning.write_bytes(b">empty\n")
    result = subprocess.run([binary, "validate", str(warning), "--json", "--error-format=json"], capture_output=True, check=False)
    assert result.returncode == 2 and result.stderr == b""
    report = json.loads(result.stdout)
    assert report["kind"] == "summary" and report["failure"] is None
    assert report["validation"]["outcome"] == "warnings"
    assert report["validation"]["codes"][0]["code"] == "fasta.empty_record"
print("diagnostic consumer: passed")
PY
