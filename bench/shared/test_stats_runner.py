#!/usr/bin/env python3
"""Check stats FAI views without launching measurements."""

import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class StatsRunnerTests(unittest.TestCase):
    def test_fai_view_preserves_shared_files_and_cleans_up(self):
        source = (ROOT / "bench/stats/run.sh").read_text()
        function = source.split("# A private FAI-only view", 1)[1].split("# Timed tools for one FASTA.", 1)[0]
        function = function[function.index("run_zebrac_tool_fai()") :]
        for mode, expected in (("success", 0), ("failure", 7), ("term", 143), ("int", 130)):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                base = Path(directory)
                fasta = base / "source with spaces.fa"
                files = {
                    fasta: b">seq\nACGT\n",
                    Path(str(fasta) + ".fai"): b"seq\t4\t5\t4\t5\n",
                    Path(str(fasta) + ".zfi"): b"original index",
                    Path(str(fasta) + ".zfi.stash"): b"prior stash",
                }
                for path, data in files.items():
                    path.write_bytes(data)
                script = r'''
set -euo pipefail
source "$1/bench/shared/runner_common.sh"
ZFASTA="$1/zig-out/bin/z-fasta"
source_fa="$2"
mode="$3"
export TMPDIR="$4"
run_zebrac_tool() {
    [[ -f "$source_fa.zfi" && -f "$source_fa.zfi.stash" ]]
    [[ -L "$view" && -L "$view.fai" && ! -e "$view.zfi" ]]
    printf '%s' "$view_dir" > "$TMPDIR/view-path"
    case "$mode" in
        failure) return 7 ;;
        term) kill -TERM "$BASHPID" ;;
        int) kill -INT "$BASHPID" ;;
        success)
            output="$(bash -c "$6")"
            [[ -z "$output" ]]
            ;;
    esac
}
'''
                script += function + '\nrun_zebrac_tool_fai "$source_fa" full tiny z-fasta-fai z-fasta out 10\n'
                result = subprocess.run(
                    ["bash", "-c", script, "bash", str(ROOT), str(fasta), mode, directory],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, expected, result.stderr)
                for path, data in files.items():
                    self.assertEqual(path.read_bytes(), data)
                view_dir = Path((base / "view-path").read_text())
                self.assertFalse(view_dir.exists())


if __name__ == "__main__":
    unittest.main()
