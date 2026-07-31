#!/usr/bin/env python3
"""Generate messy FASTA layouts into bench/shared/cache/.

Modes:
  --fixtures  tiny correctness FASTAs -> cache/messy_fixtures/
  --perf      proteome-derived layouts -> cache/messy_perf/

Default: both. Idempotent unless --force. Fixtures are the byte-level source of
truth (runners and zig tests consume the cache; nothing under messy_fixtures/ is
tracked).
"""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "bench/shared/cache"
FIXTURES_DIR = CACHE / "messy_fixtures"
PERF_DIR = CACHE / "messy_perf"
PERF_STAMP = PERF_DIR / ".stamp"
PROTEOME = ROOT / "bench/shared/data/REAL_Proteome.fasta"
MESSY_JSON = ROOT / "bench/get/messy_perf.json"
TITIN = "sp|Q8WZ42|TITIN_HUMAN"
WIDTHS = (40, 60, 70, 80, 100)
PERF_VARIANTS = ("mixed_widths", "trailing_whitespace", "mixed_crlf", "all_messy")

# Exact correctness fixtures (shared vocabulary). Edit here, then --force.
FIXTURES: dict[str, bytes] = {
    "uniform": (
        b">uniform same-width body with short final line\n"
        b"AAAACCCCGGGG\n"
        b"TTTTAAAACCCC\n"
        b"GGGGTTTT\n"
    ),
    "mixed_widths": (
        b">mixed_widths internal line widths vary\n"
        b"AAAACCCCGGGG\n"
        b"TTTTAA\n"
        b"AACCCCGGGGTT\n"
        b"TT\n"
    ),
    "trailing_whitespace": (
        b">trailing_whitespace spaces and tabs after sequence bytes\n"
        b"AAAACCCC    \n"
        b"GGGGTTTT\t\n"
        b"CCCCAAAA\n"
    ),
    "blank_lines": (
        b">blank_lines empty lines inside sequence should not count as bases\n"
        b"AAAACCCC\n"
        b"\n"
        b"GGGG\n"
        b"\n"
        b"TTTTAAAA\n"
    ),
    "mixed_crlf": (
        b">mixed_crlf mixed line endings and widths\r\n"
        b"AAAACCCC\r\n"
        b"GGGGTT\n"
        b"TTAAAA\r\n"
    ),
}


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].rstrip("\r\n")
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line))
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def wrap_mixed(seq: str, rng: random.Random) -> list[str]:
    lines: list[str] = []
    i = 0
    n = len(seq)
    while i < n:
        w = rng.choice(WIDTHS)
        chunk = seq[i : i + w]
        i += len(chunk)
        lines.append(chunk)
    return lines or [""]


def wrap_fixed(seq: str, width: int = 60) -> list[str]:
    if not seq:
        return [""]
    return [seq[i : i + width] for i in range(0, len(seq), width)]


def maybe_trail(lines: list[str], rng: random.Random, p: float = 0.3) -> list[str]:
    out: list[str] = []
    for line in lines:
        if line and rng.random() < p:
            out.append(line + "  ")
        else:
            out.append(line)
    return out


def write_records(
    path: Path,
    records: list[tuple[str, str]],
    *,
    line_fn,
    newline: str = "\n",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        for header, seq in records:
            fh.write(f">{header}{newline}")
            for line in line_fn(seq):
                fh.write(f"{line}{newline}")


def titin_length(records: list[tuple[str, str]]) -> int:
    for header, seq in records:
        name = header.split()[0]
        if name == TITIN:
            return len(seq)
    raise SystemExit(f"error: missing accession {TITIN} in {PROTEOME}")


def max_span_end(cfg: dict) -> int:
    ends: list[int] = []
    for span in cfg.get("spans", {}).values():
        region = span["region"]
        suffix = region.rsplit(":", 1)[-1]
        _, end_s = suffix.split("-", 1)
        ends.append(int(end_s))
    if not ends:
        raise SystemExit("error: messy_perf.json has no spans")
    return max(ends)


def fixtures_ready() -> bool:
    for name, data in FIXTURES.items():
        path = FIXTURES_DIR / f"{name}.fasta"
        if not path.is_file() or path.read_bytes() != data:
            return False
    return True


def generate_fixtures(*, force: bool) -> None:
    if fixtures_ready() and not force:
        print(f"messy fixtures cache ok: {FIXTURES_DIR}")
        return
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    for name, data in FIXTURES.items():
        (FIXTURES_DIR / f"{name}.fasta").write_bytes(data)
    print(f"wrote {len(FIXTURES)} messy fixtures -> {FIXTURES_DIR}")


def perf_stamp(titin_len: int, max_span_end: int) -> str:
    stat = PROTEOME.stat()
    return (
        f"proteome={PROTEOME.resolve()}\n"
        f"proteome_size={stat.st_size}\n"
        f"proteome_mtime_ns={stat.st_mtime_ns}\n"
        f"titin_len={titin_len}\n"
        f"max_span_end={max_span_end}\n"
    )


def perf_ready(titin_len: int, max_span_end: int) -> bool:
    if not PERF_STAMP.is_file():
        return False
    return (
        PERF_STAMP.read_text(encoding="utf-8") == perf_stamp(titin_len, max_span_end)
        and all((PERF_DIR / f"{name}.fasta").is_file() for name in PERF_VARIANTS)
    )


def generate_perf(*, force: bool) -> None:
    if not PROTEOME.is_file():
        raise SystemExit(
            f"error: missing {PROTEOME}; run bash bench/shared/download_data.sh"
        )
    cfg = json.loads(MESSY_JSON.read_text(encoding="utf-8"))
    need = max_span_end(cfg)
    records = parse_fasta(PROTEOME)
    tlen = titin_length(records)
    if tlen < need:
        raise SystemExit(
            f"error: {TITIN} length {tlen} < largest span end {need}; "
            "refresh bench/get/messy_perf.json spans"
        )

    if perf_ready(tlen, need) and not force:
        print(f"messy perf cache ok: {PERF_DIR}")
        return

    PERF_DIR.mkdir(parents=True, exist_ok=True)
    mixed_rng = random.Random(42)
    trail_rng = random.Random(42)
    all_rng = random.Random(43)
    trail_all_rng = random.Random(44)

    write_records(
        PERF_DIR / "mixed_widths.fasta",
        records,
        line_fn=lambda s: wrap_mixed(s, mixed_rng),
    )
    write_records(
        PERF_DIR / "trailing_whitespace.fasta",
        records,
        line_fn=lambda s: maybe_trail(wrap_fixed(s), trail_rng),
    )
    write_records(
        PERF_DIR / "mixed_crlf.fasta",
        records,
        line_fn=lambda s: wrap_fixed(s),
        newline="\r\n",
    )
    write_records(
        PERF_DIR / "all_messy.fasta",
        records,
        line_fn=lambda s: maybe_trail(wrap_mixed(s, all_rng), trail_all_rng),
    )

    PERF_STAMP.write_text(perf_stamp(tlen, need), encoding="utf-8")
    print(f"wrote {len(PERF_VARIANTS)} messy perf FASTAs -> {PERF_DIR}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--fixtures",
        action="store_true",
        help="write tiny correctness FASTAs to cache/messy_fixtures/",
    )
    ap.add_argument(
        "--perf",
        action="store_true",
        help="write proteome-derived layouts to cache/messy_perf/",
    )
    ap.add_argument("--force", action="store_true", help="rebuild even if cache exists")
    args = ap.parse_args()
    # Neither flag => both. One flag => that mode only.
    if args.fixtures or args.perf:
        do_fixtures, do_perf = args.fixtures, args.perf
    else:
        do_fixtures, do_perf = True, True
    if do_fixtures:
        generate_fixtures(force=args.force)
    if do_perf:
        generate_perf(force=args.force)


if __name__ == "__main__":
    main()
