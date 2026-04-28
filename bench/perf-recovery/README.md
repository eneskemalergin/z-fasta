# Performance Recovery Benchmarks

Focused microbenchmarks for the v0.2.6 performance recovery work.

## Startup Floor

```bash
./zig-0.16.0/zig build -Doptimize=ReleaseFast
bash bench/perf-recovery/run_startup.sh --runs 30 --warmup 5
```

The startup harness compares:

- `/bin/true`
- tiny Zig 0.16 probes built under `/tmp`
- `z-fasta --version`
- `z-fasta --help`

Results are written to `bench/perf-recovery/results/startup_<timestamp>.json`, which is ignored by git.
