# Reverse Strategy Note

This note documents the reverse-path choice currently shipped on `main` for the v0.2.8 orientation work.

## Chosen Strategy

Keep the mmap-backed reverse-streaming path as the default implementation for `z-fasta get --rc` and `--reverse-only`:

- resolve the region once using the existing O(1) indexed offset math
- start from the last base in the requested span
- walk backward through the mmapped FASTA, skipping `\n` / `\r`
- apply complement lookup on the fly when the orientation requires it
- write through the same buffered emission path used by forward extraction

This keeps the shipped reverse path zero-allocation with respect to region size and avoids creating a second copy of the extracted span.

## Why This Won For The Shipped Path

- It reuses the existing mmap-backed extraction model instead of introducing a second large-region buffering path.
- The emitted bytes stay in the same output-shaping path as forward extraction, so wrapping and header behavior remain simple to reason about.
- The measured slices show modest, predictable RC overhead relative to forward extraction without evidence that a second implementation path would pay for its complexity.
- `get --low-mem` does not currently ship, so a buffered reverse fallback would be dead surface area in v0.2.8.

## Review Harness

Focused RC review is now reproducible with `bash bench/get/run_rc_review.sh`.

- `--quick` is the default local-edit profile: smaller synthetic fixture, lighter BED batch, lighter multi-region slice.
- `--full` keeps the same benchmark shape but scales the synthetic fixture and batch sizes up for a heavier rerun.
- Outputs land in `bench/get/results/rc_review_<timestamp>/` as markdown tables, JSON, and a small RSS snapshot.

## Narrow Measurement Slice

Environment:

- host: local Linux workstation
- binary: `./zig-out/bin/z-fasta` built with `./zig build -Doptimize=ReleaseFast`
- fixture: `experiments/test_100k.fasta`, whole-sequence extraction of `seq1`
- tool: `hyperfine --warmup 3 --runs 20`

| Command | Mean [ms] | Relative |
| --- | ---: | ---: |
| `./zig-out/bin/z-fasta get experiments/test_100k.fasta seq1` | 8.2 +/- 0.7 | 1.00 |
| `./zig-out/bin/z-fasta get experiments/test_100k.fasta seq1 --rc` | 8.2 +/- 0.5 | 1.01 |
| `./zig-out/bin/z-fasta get experiments/test_100k.fasta seq1 --rc --annotate-rc` | 8.3 +/- 0.6 | 1.01 |

`/usr/bin/time -v` snapshot for the same commands:

- forward: MaxRSS `47140 kB`
- `--rc`: MaxRSS `47140 kB`
- `--rc --annotate-rc`: MaxRSS `47140 kB`

Interpretation:

- reverse-complement stayed within run-to-run noise of forward extraction on this slice
- header annotation did not produce a meaningful additional cost
- mmap dominates the observed RSS in the same way across forward and reverse runs

## Lock-In Review: Quick Profile

The quick-profile lock-in pass is not publication-grade timing, but it is fast enough for local iteration and wide enough to catch path-specific regressions.

Representative quick-profile results after the complement-validation cache fix:

| Slice | Forward | RC / related paths | Takeaway |
| --- | ---: | ---: | --- |
| large region (`250 kbp`) | `0.9 ms` | `1.1 ms` for `--rc`, `1.1 ms` for `--annotate-rc`, `0.9 ms` for `--complement-only` | RC overhead stays small on a single large region |
| full sequence (`2 Mbp`) | `3.4 ms` | `4.2 ms` for `--rc`, `4.2 ms` for `--annotate-rc` | whole-sequence RC remains in the same performance class |
| multi-region (`10` / `50`) | `0.7 ms` / `0.9 ms` | `0.8 ms` / `1.0 ms` | no-flag path stays ahead, but RC does not introduce a large multiplier |
| BED batch (`2000` rows) | `1.5 ms` | `1.9 ms` for `--rc`, `1.7 ms` for `--honor-strand --rc`, `1.5 ms` for `--complement-only` | batch complement paths are now back in the same low-ms band |

Quick-profile RSS snapshot:

- large-region forward: `1064 kB`
- large-region `--rc`: `1100 kB`
- large-region `--rc --annotate-rc`: `1084 kB`
- multi-region forward / `--rc`: `2600 kB` / `2608 kB`
- BED forward / `--honor-strand --rc`: `2112 kB` / `2104 kB`

## Optimization Finding

The review pass exposed one real overhead source worth fixing before lock-in: complement-based BED batches were repeatedly rescanning the same record sample inside the protein guard.

That is now fixed by caching the last detected record type during `ensureComplementAllowed()`. On the same quick review slice, BED `--rc` dropped from roughly `137-168 ms` down to `1.9 ms`, and BED `--complement-only` dropped from roughly `137.7 ms` down to `1.5 ms`.

## Deferred Paths

These paths remain explicitly deferred rather than shipped:

- any `get --low-mem --rc` fallback, because `get` does not currently expose a low-memory extraction mode
- a second temporary-buffer reverse path for small or medium regions, because there is not yet enough measured evidence to justify more code surface
- any extra reverse path whose only current justification would be speculation rather than a measured win

## Release Decision

For v0.2.8, keep the current mmap reverse-streaming implementation as the shipped default, keep `--annotate-rc` opt-in, keep the complement-type guard cached, and defer any extra reverse path until broader candidate benchmarking demonstrates a meaningful win.
