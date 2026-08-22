# Local benchmark tools

I want a reproducible way to re-create the benchmarks. This is still a bit of patchwork with installers, and it is currently validated on Fedora 44 x86_64. I am not pretending this is a portable package manager.

Run this from the repository root:

```bash
tools/install.sh
bash bench/shared/install_tools.sh
```

`tools/install.sh` is the one builder. When you run it, it downloads the pinned peer sources or release binaries, builds the source-based tools, compiles the Rust wrappers, installs `pyfaidx` and the report packages, checks the peer versions, strips the source-built runtime binaries, and publishes the commands under `tools/bin/`. Downloaded archives remain under `tools/build/downloads/`, extracted sources remain under `tools/src/`, and the Python environment remains under `tools/venv/`. Temporary worktrees, compiler outputs, Cargo and pip caches, and extracted release copies are removed automatically when the script exits. Set `KEEP_BUILD_ARTIFACTS=1` when you need to inspect failed build work. The benchmark scripts use `tools/bin/` by default, but you can override an individual command with variables such as `SAMTOOLS`, `BEDTOOLS`, or `SEQKIT`.

The current peer set is [samtools](https://github.com/samtools/samtools) 1.24 with [HTSlib](https://github.com/samtools/htslib) 1.24, [bedtools](https://github.com/arq5x/bedtools2) 2.31.1, [seqkit](https://github.com/shenwei356/seqkit) 2.13.0, [seqtk](https://github.com/lh3/seqtk) 1.5-r133, [fastahack](https://github.com/ekg/fastahack) 1.0.0, and [pyfaidx](https://github.com/mdshw5/pyfaidx) 0.9.0.4.

The Rust benchmark lanes are small CLI adapters, not renamed upstream programs. The `noodles` command uses [noodles](https://github.com/zaeleus/noodles) and reports noodles-fasta 0.66.0. The `rustbio` command uses [rust-bio](https://github.com/rust-bio/rust-bio) 4.0.1 and keeps its custom strict FAI indexer. Both expose `index`, `get`, and `stats` so the benchmark can exercise the same basic operations. Their `stats` subcommands use the shared `tools/stats_peer.rs` formulas; that file is compiled into the wrappers and is not a separate runtime tool.

The benchmark runner is [zebrac](https://github.com/eneskemalergin/zebrac), kept at `tools/zebrac`. It is separate from the peer installer and is not built by `tools/install.sh`. The local binary currently reports zebrac 0.6.2.

The second command only verifies that the local bundle is complete and that the commands resolve to the expected versions. It does not silently fall back to tools installed elsewhere on the system.
