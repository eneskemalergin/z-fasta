#!/usr/bin/env python3
"""Generate synthetic scaling FASTAs into bench/shared/cache/scaling/.

One shared tree for index and stats (index uses size + budget + fixed; stats uses
size + fixed only). Generate-if-missing unless --force. Stamp fingerprints params.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "bench/shared/cache/scaling"
STAMP = OUT_DIR / ".stamp"

LINE = "ACGTACGTAC" * 8
BUDGET_BYTES = 50 * 1024 * 1024
MIN_SEQ_LEN = 80
FIXED_SEQ_LEN = 1024
SIZE_MBS = (1, 5, 10, 25, 50, 100, 250, 500, 1000)
BUDGET_COUNTS = (1_000, 10_000, 100_000, 250_000)
FIXED_COUNTS = (100_000, 250_000, 500_000, 1_000_000)
SIZE_SEQ_COUNT = 100
SCHEMA = "scaling.v1"


def write_wrapped(handle, seq_len: int) -> None:
    remaining = seq_len
    while remaining > 0:
        chunk = min(80, remaining)
        handle.write(LINE[:chunk] + "\n")
        remaining -= chunk


def write_budget(count: int, path: Path) -> None:
    seq_len = max(MIN_SEQ_LEN, BUDGET_BYTES // count - 20)
    with path.open("w", encoding="ascii") as handle:
        for i in range(1, count + 1):
            handle.write(f">seq{i}\n")
            write_wrapped(handle, seq_len)


def write_fixed(count: int, path: Path) -> None:
    with path.open("w", encoding="ascii") as handle:
        for i in range(1, count + 1):
            handle.write(f">seq{i}\n")
            write_wrapped(handle, FIXED_SEQ_LEN)


def write_size(mb: int, path: Path) -> None:
    total = mb * 1024 * 1024
    with path.open("w", encoding="ascii") as handle:
        for i in range(1, SIZE_SEQ_COUNT + 1):
            handle.write(f">seq{i}\n")
            # -50 matches legacy inline generator (100 seqs, ~target MiB on disk)
            write_wrapped(handle, max(MIN_SEQ_LEN, total // SIZE_SEQ_COUNT - 50))


def expected_paths(mode: str) -> list[Path]:
    paths: list[Path] = []
    if mode in ("all", "size"):
        paths.extend(OUT_DIR / f"size_{mb}mb.fasta" for mb in SIZE_MBS)
    if mode in ("all", "seq", "budget"):
        paths.extend(OUT_DIR / f"seqs_budget_{c}.fasta" for c in BUDGET_COUNTS)
    if mode in ("all", "seq", "fixed"):
        paths.extend(OUT_DIR / f"seqs_fixed_{c}.fasta" for c in FIXED_COUNTS)
    return paths


def stamp_payload() -> str:
    parts = [
        f"schema={SCHEMA}",
        f"size_mbs={','.join(map(str, SIZE_MBS))}",
        f"budget_counts={','.join(map(str, BUDGET_COUNTS))}",
        f"fixed_counts={','.join(map(str, FIXED_COUNTS))}",
        f"size_seq_count={SIZE_SEQ_COUNT}",
        f"budget_bytes={BUDGET_BYTES}",
        f"fixed_seq_len={FIXED_SEQ_LEN}",
        f"min_seq_len={MIN_SEQ_LEN}",
        f"line_hash={hashlib.sha256(LINE.encode()).hexdigest()[:16]}",
    ]
    return "\n".join(parts) + "\n"


def stamp_ok() -> bool:
    if not STAMP.is_file():
        return False
    return STAMP.read_text(encoding="utf-8") == stamp_payload()


def outputs_ready(mode: str) -> bool:
    if not stamp_ok():
        return False
    return all(p.is_file() for p in expected_paths(mode))


def clean_legacy() -> None:
    """Remove pre-budget seqs_*.fasta names that are neither budget nor fixed."""
    if not OUT_DIR.is_dir():
        return
    for path in OUT_DIR.glob("seqs_*.fasta"):
        if path.name.startswith(("seqs_budget_", "seqs_fixed_")):
            continue
        path.unlink()


def generate(*, mode: str, force: bool, do_clean_legacy: bool) -> None:
    if mode not in ("all", "size", "seq", "budget", "fixed"):
        raise SystemExit(f"error: mode must be all|size|seq|budget|fixed, got {mode!r}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if do_clean_legacy:
        clean_legacy()

    if outputs_ready(mode) and not force:
        print(f"scaling cache ok: {OUT_DIR} ({mode})")
        return

    if mode in ("all", "size"):
        for mb in SIZE_MBS:
            path = OUT_DIR / f"size_{mb}mb.fasta"
            if force or not path.is_file():
                write_size(mb, path)
    if mode in ("all", "seq", "budget"):
        for count in BUDGET_COUNTS:
            path = OUT_DIR / f"seqs_budget_{count}.fasta"
            if force or not path.is_file():
                write_budget(count, path)
    if mode in ("all", "seq", "fixed"):
        for count in FIXED_COUNTS:
            path = OUT_DIR / f"seqs_fixed_{count}.fasta"
            if force or not path.is_file():
                write_fixed(count, path)

    # Stamp always reflects full param set (not mode subset).
    STAMP.write_text(stamp_payload(), encoding="utf-8")
    print(f"wrote scaling FASTAs -> {OUT_DIR} ({mode})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--mode",
        choices=("all", "size", "seq", "budget", "fixed"),
        default="all",
        help="which fixture families to ensure (seq=budget+fixed; default: all)",
    )
    ap.add_argument("--force", action="store_true", help="rebuild even if present")
    ap.add_argument(
        "--clean-legacy",
        action="store_true",
        help="delete obsolete seqs_*.fasta that are not budget/fixed",
    )
    args = ap.parse_args()
    generate(mode=args.mode, force=args.force, do_clean_legacy=args.clean_legacy)


if __name__ == "__main__":
    main()
