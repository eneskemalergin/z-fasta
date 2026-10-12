#!/usr/bin/env python3
"""Check benchmark readers with small JSON fixtures, without measurements."""

import importlib.util
import json
import subprocess
import os
import copy
import sys
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def load_report(suite):
    spec = importlib.util.spec_from_file_location(
        f"report_{suite}", ROOT / "bench" / suite / "generate_report.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReportReaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = load_report("index")
        cls.get = load_report("get")
        cls.stats = load_report("stats")

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="zfasta-report-test-")
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "sample.json"

    def result(self):
        return {
            "command": "bash -c 'echo hi'",
            "argv": ["bash", "-c", "echo hi"],
            "wall_time": {"mean": 1_000_000_000},
            "peak_rss": {"mean": 1_048_576},
        }

    def read(self, suite, result, metadata=None):
        self.path.write_text(json.dumps({"zebrac_version": "0.6.2", "results": [result]}))
        if suite == "index":
            return self.index.load_zebrac_json(self.path, metadata)
        return getattr(self, suite).load_zebrac_json(self.path, metadata, self.index)

    def test_all_readers_reject_undeclared_result_formats(self):
        self.path.write_text(json.dumps({"results": [self.result()]}))
        for suite in ("index", "get", "stats"):
            with self.subTest(suite=suite):
                with self.assertRaisesRegex(ValueError, "unsupported format"):
                    if suite == "index":
                        self.index.load_zebrac_json(self.path)
                    else:
                        getattr(self, suite).load_zebrac_json(self.path, None, self.index)

    def test_all_readers_reject_malformed_or_empty_result_lists(self):
        invalid = [{"zebrac_version": "0.6.2"}]
        for results in (None, {}, "", [], [None], [True], [1], ["command"], [[]]):
            invalid.append({"zebrac_version": "0.6.2", "results": results})
        for suite in ("index", "get", "stats"):
            for data in invalid:
                with self.subTest(suite=suite, data=data):
                    self.path.write_text(json.dumps(data))
                    with self.assertRaisesRegex(ValueError, "results"):
                        if suite == "index":
                            self.index.load_zebrac_json(self.path)
                        else:
                            getattr(self, suite).load_zebrac_json(self.path, None, self.index)

    def test_valid_means_keep_seconds_and_mib(self):
        for suite in ("index", "get", "stats"):
            with self.subTest(suite=suite):
                row = self.read(suite, self.result()).iloc[0]
                self.assertEqual(row["mean"], 1)
                self.assertEqual(row["peak_rss_mb"], 1)

    def test_failed_samples_never_enter_report_tables(self):
        for suite in ("index", "get", "stats"):
            for value in (1, -1, True, "0", None, 0.5):
                with self.subTest(suite=suite, failed_sample_count=value):
                    result = self.result()
                    result["failed_sample_count"] = value
                    with self.assertRaisesRegex(ValueError, "failed_sample_count"):
                        self.read(suite, result)
            for value in (0, -1, True, "1", None, 0.5):
                with self.subTest(suite=suite, sample_count=value):
                    result = self.result()
                    result["sample_count"] = value
                    with self.assertRaisesRegex(ValueError, "sample_count"):
                        self.read(suite, result)
            result = self.result()
            result.update(sample_count=1, failed_sample_count=0)
            self.assertEqual(self.read(suite, result).iloc[0]["mean"], 1)

    def test_declared_metric_units_match_report_conversions(self):
        for suite in ("index", "get", "stats"):
            for metric, expected, wrong in (("wall_time", "nanoseconds", "seconds"),
                                            ("peak_rss", "bytes", "kilobytes")):
                for unit in (wrong, None, "", 1):
                    with self.subTest(suite=suite, metric=metric, unit=unit):
                        result = self.result()
                        result[metric]["unit"] = unit
                        with self.assertRaisesRegex(ValueError, metric + r"\.unit"):
                            self.read(suite, result)
                result = self.result()
                result[metric]["unit"] = expected
                row = self.read(suite, result).iloc[0]
                self.assertEqual(row["mean"], 1)
                self.assertEqual(row["peak_rss_mb"], 1)

    def test_missing_or_invalid_means_never_become_zero(self):
        for suite in ("index", "get", "stats"):
            for metric in ("wall_time", "peak_rss"):
                for value in (None, False, -1, float("nan"), float("inf"), "0"):
                    with self.subTest(suite=suite, metric=metric, value=value):
                        result = self.result()
                        result[metric] = {"mean": value} if value is not None else {}
                        with self.assertRaisesRegex(ValueError, metric + r"\.mean"):
                            self.read(suite, result)
                result = self.result()
                del result[metric]
                with self.assertRaisesRegex(ValueError, metric + r"\.mean"):
                    self.read(suite, result)

    def metadata(self, *commands):
        return self.index.pd.DataFrame([
            {"raw_json": str(self.path), "command": command, "tool": tool}
            for command, tool in commands
        ])

    def test_argv_matches_reordered_metadata_without_position_fallback(self):
        for suite in ("index", "get", "stats"):
            result = self.result()
            result["command"] = r"bash -c echo\ hi"
            metadata = self.metadata(("bash -c 'echo other'", "wrong"), ("bash -c 'echo hi'", "wanted"))
            with self.subTest(suite=suite):
                self.assertEqual(self.read(suite, result, metadata).iloc[0]["tool"], "wanted")

    def test_unmatched_or_conflicting_metadata_is_rejected(self):
        cases = [
            self.metadata(("bash -c 'echo other'", "wrong")),
            self.metadata(("bash -c 'echo hi'", "first"), ("bash -c 'echo hi'", "second")),
        ]
        for suite in ("index", "get", "stats"):
            for metadata in cases:
                with self.subTest(suite=suite, rows=len(metadata)):
                    with self.assertRaisesRegex(ValueError, "metadata match"):
                        self.read(suite, self.result(), metadata)

    def test_supplied_metadata_requires_rows_for_each_result_file(self):
        other_file = self.metadata(("bash -c 'echo hi'", "wrong"))
        other_file["raw_json"] = str(self.path.parent / "other.json")
        missing_path = self.metadata(("bash -c 'echo hi'", "wrong")).drop(columns="raw_json")
        for suite in ("index", "get", "stats"):
            for metadata in (other_file, missing_path, self.index.pd.DataFrame()):
                with self.subTest(suite=suite, columns=list(metadata.columns), rows=len(metadata)):
                    with self.assertRaisesRegex(ValueError, "metadata match"):
                        self.read(suite, self.result(), metadata)
            with self.subTest(suite=suite, metadata="absent legacy metadata"):
                self.assertEqual(self.read(suite, self.result()).iloc[0]["tool"],
                                 {"index": self.result()["command"], "get": "", "stats": "unknown"}[suite])

    def test_declared_metadata_is_not_treated_as_absent(self):
        directory = Path(self.directory.name)
        manifest = {"metadata": "metadata_test.jsonl"}
        with self.assertRaisesRegex(ValueError, "metadata file"):
            self.index.load_metadata(directory, manifest)
        (directory / manifest["metadata"]).write_text(" \n")
        metadata = self.index.load_metadata(directory, manifest)
        self.assertIsNotNone(metadata)
        self.assertTrue(metadata.empty)
        for suite in ("index", "get", "stats"):
            with self.subTest(suite=suite):
                with self.assertRaisesRegex(ValueError, "metadata match"):
                    self.read(suite, self.result(), metadata)
        self.assertIsNone(self.index.load_metadata(directory, {}))
        (directory / manifest["metadata"]).write_text(json.dumps({
            "raw_json": str(self.path), "command": self.result()["command"], "tool": "wanted"
        }) + "\n")
        metadata = self.index.load_metadata(directory, manifest)
        self.assertEqual(self.read("index", self.result(), metadata).iloc[0]["tool"], "wanted")

    def test_malformed_metadata_declarations_are_not_absent(self):
        directory = Path(self.directory.name)
        for declaration in (None, False, 0, [], {}, ""):
            with self.subTest(metadata=declaration):
                with self.assertRaisesRegex(ValueError, "metadata"):
                    self.index.load_metadata(directory, {"metadata": declaration})
        self.assertIsNone(self.index.load_metadata(directory, {}))

    def test_latest_never_borrows_an_older_manifest_or_metadata_file(self):
        directory = Path(self.directory.name)
        old = directory / "run_20260101.json"
        old.write_text(json.dumps({"timestamp": "20260101"}))
        old_metadata = directory / "metadata_20260101.jsonl"
        old_metadata.write_text(json.dumps({"tool": "old"}) + "\n")
        latest = directory / "LATEST"
        latest.write_text("20260201\n")
        with self.assertRaisesRegex(ValueError, "LATEST"):
            self.index.load_run_manifest(directory)
        self.assertIsNone(self.index.load_metadata(directory))
        latest.write_text("\n")
        with self.assertRaisesRegex(ValueError, "LATEST"):
            self.index.load_run_manifest(directory)
        with self.assertRaisesRegex(ValueError, "LATEST"):
            self.index.load_metadata(directory)
        latest.write_text("20260201\n")
        current = directory / "run_20260201.json"
        current.write_text(json.dumps({"timestamp": "20260201"}))
        (directory / "metadata_20260201.jsonl").write_text(json.dumps({"tool": "current"}) + "\n")
        self.assertEqual(self.index.load_run_manifest(directory)["timestamp"], "20260201")
        self.assertEqual(self.index.load_metadata(directory).iloc[0]["tool"], "current")
        self.assertEqual(self.index.load_run_manifest(directory, old)["timestamp"], "20260101")
        latest.unlink()
        self.assertEqual(self.index.load_run_manifest(directory)["timestamp"], "20260201")
        self.assertEqual(self.index.load_metadata(directory).iloc[0]["tool"], "current")

    def test_declared_sections_must_exist_as_the_expected_directory(self):
        directory = Path(self.directory.name)
        section_file = directory / "not-a-directory"
        section_file.write_text("")
        valid = directory / "valid"
        valid.mkdir()
        (valid / "Genome.json").write_text(json.dumps({
            "zebrac_version": "0.6.2", "results": [self.result()]
        }))
        for suite, key in (("index", "real"), ("get", "perf_pos"), ("stats", "perf_full")):
            def load(manifest):
                if suite == "index":
                    return self.index.load_perf(directory, None, manifest)
                return getattr(self, suite).load_section(directory, manifest, key, self.index)

            for name in ("missing", section_file.name):
                with self.subTest(suite=suite, section=name):
                    with self.assertRaisesRegex(ValueError, "declared.*section"):
                        load({"sections": {key: name}})
            with self.subTest(suite=suite, section="present"):
                self.assertEqual(len(load({"sections": {key: valid.name}})), 1)
            with self.subTest(suite=suite, section="omitted"):
                self.assertIsNone(load({"sections": {}}))

    def test_identical_metadata_duplicates_keep_the_same_identity(self):
        metadata = self.metadata(("bash -c 'echo hi'", "wanted"), ("bash -c 'echo hi'", "wanted"))
        metadata["output_bases"] = [None, None]
        extra = self.metadata(("bash -c 'echo other'", "other"))
        extra["output_bases"] = 100
        metadata = self.index.pd.concat([metadata, extra], ignore_index=True)
        for suite in ("index", "get", "stats"):
            with self.subTest(suite=suite):
                self.assertEqual(self.read(suite, self.result(), metadata).iloc[0]["tool"], "wanted")

    def test_get_throughput_uses_million_bases_per_second(self):
        metadata = self.metadata(("bash -c 'echo hi'", "wanted"))
        metadata["output_bases"] = 1_000_000
        row = self.read("get", self.result(), metadata).iloc[0]
        self.assertEqual(row["throughput_mbp_s"], 1)
        self.assertNotIn("throughput_mibs", row.index)

    def test_build_facts_survive_reading_and_report_text(self):
        identity = {"checkout_revision": "first", "checkout_dirty": False,
                    "executable": {"size": 10}, "declared_build": {"optimize": None}}
        metadata = self.metadata(("bash -c 'echo hi'", "wanted"))
        metadata["z_fasta_build"] = [identity]
        manifest = {"z_fasta_build": identity}
        for suite in ("index", "get", "stats"):
            with self.subTest(suite=suite):
                self.assertEqual(self.read(suite, self.result(), metadata).iloc[0]["z_fasta_build"], identity)
                report = getattr(self, suite)
                text = report.md_run_provenance(manifest) if suite == "index" else report.md_run_provenance(manifest, self.index)
                self.assertIn('"checkout_revision": "first"', text)
        self.assertEqual(self.index.md_build_identity({"z_fasta": "0.3.4"}), "")

    def test_build_facts_show_no_absolute_executable_path(self):
        for path, shown in ((ROOT / "zig-out/bin/z-fasta", "zig-out/bin/z-fasta"), (Path("/elsewhere/bin/z-fasta"), "z-fasta")):
            with self.subTest(path=path):
                identity = {"checkout_revision": "first", "executable": {"path": str(path), "size": 10}}
                text = self.index.md_build_identity({"z_fasta_build": identity})
                self.assertIn(f'"path": "{shown}"', text)
                self.assertIn('"size": 10', text)
                self.assertEqual(identity["executable"]["path"], str(path))

    def test_stats_resume_requires_unchanged_recorded_files(self):
        source = (ROOT / "bench/stats/run.sh").read_text()
        function = source.split("write_run_manifest() {", 1)[1]
        code = function.split("<<'PY'\n", 1)[1].split("\nPY\n}", 1)[0]
        identity = {"subject": {"checkout_revision": "first", "checkout_dirty": False,
                    "executable": {"size": 10, "mtime_ns": 1},
                    "declared_build": {"optimize": "ReleaseFast"}},
                    "inputs": {"fixture": {"size": 4}}, "tools": {"samtools": {"size": 10}}}
        cases = [("same", identity, True)]
        for field in ("revision", "dirty", "executable", "build", "input", "tool"):
            changed = copy.deepcopy(identity)
            if field == "revision":
                changed["subject"]["checkout_revision"] = "second"
            elif field == "dirty":
                changed["subject"]["checkout_dirty"] = True
            elif field == "executable":
                changed["subject"]["executable"]["mtime_ns"] = 2
            elif field == "build":
                changed["subject"]["declared_build"]["optimize"] = "Debug"
            elif field == "input":
                changed["inputs"]["fixture"]["size"] = 5
            else:
                changed["tools"]["samtools"]["size"] = 11
            cases.append((field, changed, False))
        for field in ("unknown_tool", "no_tools", "missing_tools"):
            unknown = copy.deepcopy(identity)
            if field == "unknown_tool":
                unknown["tools"]["samtools"] = None
            elif field == "no_tools":
                unknown["tools"] = {}
            else:
                del unknown["tools"]
            cases.append((field, unknown, False))
        manifest = Path(self.directory.name) / "run_test.json"
        (manifest.parent / "verify_test.log").write_text("Results: 57 passed, 0 failed\nALL PASSED\n")
        for name, current, expected in cases:
            with self.subTest(change=name):
                prior_identity = current if name in ("unknown_tool", "no_tools", "missing_tools") else identity
                manifest.write_text(json.dumps({"verification_identity": prior_identity,
                                    "verify_skipped": False, "verify_pass": "57"}))
                env = dict(os.environ, BENCH_VERIFICATION_MODULE_DIR=str(ROOT / "bench/shared"),
                           BENCH_ZFASTA_BUILD=json.dumps(current["subject"]),
                           BENCH_CURRENT_VERIFICATION=json.dumps(current))
                subprocess.run([sys.executable, "-c", code, str(manifest), "test", "meta", "1", "0", "1",
                                "false", "full", "", ""], env=env, check=True, capture_output=True)
                result = json.loads(manifest.read_text())
                self.assertEqual(result["verify_skipped"], not expected)
                self.assertEqual(self.stats.resolve_verify_status(result, manifest.parent),
                                 (False, "57") if expected else (True, None))
        manifest.write_text(json.dumps({"verify_skipped": False, "verify_pass": "57"}))
        env["BENCH_CURRENT_VERIFICATION"] = json.dumps(identity)
        subprocess.run([sys.executable, "-c", code, str(manifest), "test", "meta", "1", "0", "1",
                        "false", "full", "", ""], env=env, check=True, capture_output=True)
        self.assertEqual(self.stats.resolve_verify_status(json.loads(manifest.read_text()), manifest.parent), (True, None))

    def region_bases(self, region, fasta=""):
        source = (ROOT / "bench/get/run.sh").read_text()
        function = source.split("region_span_bases() {", 1)[1]
        code = function.split("<<'PY'\n", 1)[1].split("\nPY\n}", 1)[0]
        return subprocess.run([sys.executable, "-c", code, region, str(fasta)], capture_output=True, text=True)

    def test_region_base_counts_follow_indexed_bounds_and_literal_names(self):
        fasta = Path(self.directory.name) / "regions.fa"
        data = bytearray()
        rows = []
        for name, sequence in (("seq", "ACGT"), ("long", "A" * 16), ("name:1-2", "A" * 6), ("seq:18446744073709551616-1", "AA"), ("seq:1_-2", "AAA")):
            header = f">{name}\n".encode()
            offset = len(data) + len(header)
            data.extend(header + sequence.encode() + b"\n")
            rows.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n")
        fasta.write_bytes(data)
        Path(str(fasta) + ".fai").write_text("".join(rows))
        cases = {"seq": 4, "seq:1-999": 4, "seq:+1-+3": 3, "long:1__0-1__2": 3, "name:1-2:1-3": 3, "seq:18446744073709551616-1": 2, "seq:1_-2": 3, "seq:1-": 4}
        for region, expected in cases.items():
            with self.subTest(region=region):
                result = self.region_bases(region, fasta)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(int(result.stdout), expected)
        for region in ("seq:0-1", "seq:5-", "seq:3-2", "seq:1--0", " seq"):
            with self.subTest(region=region):
                self.assertNotEqual(self.region_bases(region, fasta).returncode, 0)

    def test_unknown_open_end_never_invents_one_output_base(self):
        self.assertNotEqual(self.region_bases("seq:1-").returncode, 0)
        explicit = self.region_bases("seq:1-3")
        self.assertEqual(explicit.returncode, 0, explicit.stderr)
        self.assertEqual(int(explicit.stdout), 3)


if __name__ == "__main__":
    unittest.main()
