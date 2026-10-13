#!/usr/bin/env python3
"""Generate bench/validate/REPORT.md from the latest run's correctness table and zebrac JSON.

Timing rows are read through the stats report's zebrac loader, which builds on the index
report's reader, so every suite reads zebrac output the same way. Tables only, no figures.
"""

from __future__ import annotations

import csv
import importlib.util
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
REPORT_PATH = SCRIPT_DIR / "REPORT.md"
DATASETS = ["Genome", "Transcriptome", "Proteome"]
TOOLS = ["z-fasta", "z-fasta-index", "seqkit", "samtools"]
TOOL_LABELS = {
    "z-fasta": "z-fasta validate",
    "z-fasta-index": "z-fasta validate --index",
    "seqkit": "SeqKit seq -v (partial)",
    "samtools": "samtools faidx (partial)",
}
GROUP_LABELS = {
    "crafted": "Crafted proteome copies",
    "fixture": "Messy fixtures",
    "variant": "Messy proteome variants",
    "real": "REAL datasets",
}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_tests(manifest: dict) -> list[dict]:
    name = manifest.get("sections", {}).get("correctness")
    if not name:
        return []
    path = RESULTS_DIR / name
    if not path.is_file():
        raise ValueError(f"missing declared correctness section: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def describe(value: str) -> str:
    """Spells out an `exit N; code,code` verdict for a list item."""
    match = re.fullmatch(r"exit (\d+); (.*)", value)
    if not match:
        return value
    codes = match.group(2)
    found = "no findings" if codes == "none" else ", ".join(f"`{code}`" for code in codes.split(","))
    return f"exit {match.group(1)} with {found}"


def md_overview(manifest: dict, tests: list[dict]) -> str:
    passed = sum(row["status"] == "pass" for row in tests)
    failed = sum(row["status"] == "fail" for row in tests)
    if tests:
        verdict = f"Correctness reported **{passed}** passing and **{failed}** failing checks before performance (see `bench/validate/run.sh`)."
    else:
        verdict = "Correctness was skipped for this run; run the suite without `--skip-tests` before relying on these numbers."
    return "\n".join([
        "## Overview",
        "",
        "This report checks `z-fasta validate` on crafted, messy, and REAL inputs, then times it beside two partial peers on the shared REAL datasets (Genome, Transcriptome, Proteome).",
        "",
        "**What each tool checks**",
        "",
        "- **z-fasta validate:** structure, identifiers, alphabet, and line layout in one streaming pass through a 256 KiB buffer, with every occurrence counted. Timed as `validate --json`, which prints one report object; the alphabet is chosen automatically from the first 100,000 symbols.",
        "- **z-fasta validate --index:** the same report, plus a record-by-record comparison with the `.zfi` built beforehand, made in the same pass by `index`'s scanner.",
        "- **SeqKit `seq -v`** (partial): checks every sequence against an explicit alphabet (`-t dna` or `-t protein`) on one thread (`-j 1`) and writes the records back out (`-o /dev/null`), so its time includes that rewrite. It does not check layout or headers.",
        "- **samtools `faidx`** (partial): builds a `.fai`, which needs fixed line geometry. It checks layout, not the alphabet.",
        "",
        verdict,
    ])


def md_correctness(tests: list[dict]) -> str:
    if not tests:
        return ""
    checked = [row for row in tests if row["status"] != "info"]
    groups: dict[tuple[str, str], list[int]] = {}
    for row in checked:
        key = (row["group"], row["check"])
        counts = groups.setdefault(key, [0, 0])
        counts[0] += row["status"] == "pass"
        counts[1] += 1
    summary = table(
        ["Inputs", "Check", "Passed"],
        [[GROUP_LABELS.get(group, group), check, f"{ok} of {total}"] for (group, check), (ok, total) in groups.items()],
    )

    by_case: dict[tuple[str, str], dict[str, dict]] = {}
    for row in checked:
        by_case.setdefault((row["group"], row["case"]), {})[row["check"]] = row
    verdicts: dict[str, list[str]] = {}
    for (group, case), checks in by_case.items():
        if group not in ("real", "fixture", "variant"):
            continue
        findings = checks.get("findings", {}).get("observed", "")
        accepted = checks.get("accepted (index, validate)", {}).get("observed", "")
        parts = [f"validate {describe(findings)}" if findings else (f"validate {accepted}, findings not asserted" if accepted else "validate not checked")]
        for check, label in (("accepted (index, validate)", "index"), ("index pair (validate --index)", "validate --index"), ("layout (samtools, validate)", "samtools faidx"), ("alphabet (seqkit, validate)", "SeqKit seq -v")):
            row = checks.get(check)
            if row:
                verdict = row["observed"] if check == "index pair (validate --index)" else row["expected"]
                parts.append(f"{label} {verdict}")
        verdicts.setdefault(GROUP_LABELS[group], []).append(f"  - `{case}`: " + "; ".join(parts) + ".")
    crafted = []
    for row in checked:
        if row["group"] != "crafted" or row["check"] != "findings":
            continue
        if row["status"] == "pass":
            crafted.append(f"- `{row['case']}`: {describe(row['observed'])}, as expected.")
        else:
            crafted.append(f"- `{row['case']}`: expected {describe(row['expected'])}, observed {describe(row['observed'])}; **failed**.")
    verdict_lines = []
    for label, items in verdicts.items():
        verdict_lines.append(f"- {label}:")
        verdict_lines += items
    differences = [
        [f"`{row['case']}`", row["expected"], row["observed"]]
        for row in tests
        if row["status"] == "info"
    ]
    parts = [
        "## Correctness",
        "",
        "Expected findings come from how each input is built and from the documented exit rules: 0 for a clean file, 2 for warnings only, and 1 for errors. A file counts as accepted when validate reports no error-level code.",
        "",
        summary,
        "",
        "### Proteome copies with one injected defect",
        "",
        "\n".join(crafted),
        "",
        "### Verdicts on real and messy files",
        "",
        "validate's findings, with the verdict each other tool gives the same file. `index` accepts every one, and validate accepts the same files; `validate --index` then matches every index built. samtools rejects exactly the files where validate reports a wrapping, whitespace, blank-line, or line-ending finding.",
        "",
        "\n".join(verdict_lines),
    ]
    if differences:
        parts += [
            "",
            "### Rules that differ (not asserted)",
            "",
            "Each file holds one byte in question inside 240 nucleotide bases, checked with SeqKit `-t dna`. These verdicts are recorded each run and not asserted, since the two alphabets are defined differently.",
            "",
            table(["Byte in question", "SeqKit seq -v", "validate"], differences),
        ]
    return "\n".join(parts)


def md_performance(perf) -> str:
    if perf is None or perf.empty:
        return "## Performance\n\nNo timing results in this run."
    rows_time, rows_rss, rows_faults, rows_rate = [], [], [], []
    for dataset in DATASETS:
        subset = perf[perf["workload"] == dataset]
        if subset.empty:
            continue
        cells = {"time": [dataset], "rss": [dataset], "faults": [dataset], "rate": [dataset]}
        for tool in TOOLS:
            row = subset[subset["tool"] == tool]
            if row.empty:
                for key in cells:
                    cells[key].append("")
                continue
            row = row.iloc[0]
            cells["time"].append(f"{row['mean']:.3f} +/- {row['stddev']:.3f} s")
            cells["rss"].append(f"{row['peak_rss_mb']:,.1f} MiB")
            cells["faults"].append(f"{row['minor_faults']:,.0f}")
            rate = row["throughput_mibs"]
            cells["rate"].append("" if rate is None or rate != rate else f"{rate:,.0f} MiB/s")
        rows_time.append(cells["time"])
        rows_rss.append(cells["rss"])
        rows_faults.append(cells["faults"])
        rows_rate.append(cells["rate"])
    header = ["Dataset"] + [TOOL_LABELS[tool] for tool in TOOLS]
    return "\n".join([
        "## Performance",
        "",
        "zebrac times each command in turn, in one batch, on a warm page cache. The peers check less than validate (see the overview), so their times are references, not a ranking.",
        "",
        "### Wall time (mean +/- standard deviation)",
        "",
        table(header, rows_time),
        "",
        "### Peak memory (mean resident set)",
        "",
        "validate reads through a fixed buffer, so its peak memory does not grow with sequence bytes; its duplicate check keeps a fingerprint and an offset per identifier. With `--index` it reads the index in file order through fixed buffers.",
        "",
        table(header, rows_rss),
        "",
        "### Minor page faults (mean)",
        "",
        table(header, rows_faults),
        "",
        "### Throughput (input MiB per second of wall time)",
        "",
        table(header, rows_rate),
    ])


def md_provenance(manifest: dict, ir) -> str:
    tools = manifest.get("tools", {})
    lines = [
        "## Run Provenance",
        "",
        f"- Run: `{manifest.get('timestamp', '')}`, zebrac `{manifest.get('zebrac', '')}`, {manifest.get('runs')} measured samples, {manifest.get('warmup')} warmups, {manifest.get('duration_ms')} ms budget per command.",
        f"- z-fasta: `{manifest.get('z_fasta', '')}`.",
        f"- samtools: `{tools.get('samtools', 'unknown')}`; SeqKit: `{tools.get('seqkit', 'unknown')}`.",
    ]
    identity = ir.md_build_identity(manifest)
    if identity:
        lines += ["", identity]
    return "\n".join(lines)


def main() -> None:
    ir = load_module("index_generate_report", SCRIPT_DIR.parent / "index" / "generate_report.py")
    sr = load_module("stats_generate_report", SCRIPT_DIR.parent / "stats" / "generate_report.py")
    manifest = ir.load_run_manifest(RESULTS_DIR)
    if manifest is None:
        raise SystemExit("error: no validate run manifest; run bash bench/validate/run.sh first")
    tests = load_tests(manifest)
    perf = sr.load_section(RESULTS_DIR, manifest, "perf", ir)
    parts = [
        "<!-- markdownlint-disable MD013 MD041 -->",
        "",
        "# z-fasta Validate Benchmark Report",
        "",
        "_Auto-generated by `bench/validate/generate_report.py` from the correctness table and zebrac results._",
        "",
        md_overview(manifest, tests),
        "",
        md_correctness(tests),
        "",
        md_performance(perf),
        "",
        md_provenance(manifest, ir),
        "",
    ]
    REPORT_PATH.write_text(ir.normalize_markdown(parts))
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
