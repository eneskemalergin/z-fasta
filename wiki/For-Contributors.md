# For Contributors

I prefer changes that start from observable behavior and then trace inward to the implementation. z-fasta shares contracts across commands, so a small parser change can reach farther than one source file.

## Build and test

Install Zig 0.16.0 from the [official Zig download page](https://ziglang.org/download/) and confirm the compiler before building:

```bash
zig version
zig build
zig build test --summary all
zig build test -Doptimize=ReleaseFast --summary all
zig build -Doptimize=ReleaseFast
```

The version command must print `0.16.0`. A compiler installed on your `PATH` is fine when it is that exact version. The repository does not provide a public Zig wrapper or compiler download.

## Source map

- `src/main.zig`: CLI parsing, help, dispatch, and diagnostics.
- `src/indexer.zig`: bounded FASTA scanner and index output.
- `src/index_format.zig`: `.zfi` or `.fai` loading, validation, lookup, and ownership.
- `src/getter.zig`: requests, coordinates, transforms, and bounded output.
- `src/stats.zig`: length and composition reports plus shared type detection.
- `src/validator.zig`: events, JSON, fix policy, and normalized rewrite.
- `src/bed_parser.zig`: BED row conversion and strand parsing.
- `src/complement.zig`: IUPAC complement table.
- `src/platform.zig`: portable read-only mapping and advice.

The Wiki owns public behavior and workflow documentation. Source files and their tests own implementation details that are not part of the public contract.

## Local correctness gates

```bash
bash bench/index/run.sh --skip-benchmarks --skip-messy --skip-report
bash bench/get/run.sh --skip-benchmarks --skip-report
bash bench/stats/run.sh --skip-benchmarks --skip-report
```

Run every affected gate after code changes. Synthetic unit fixtures do not replace real-data and peer verification.

## Contracts to preserve

- A present `.zfi` is authoritative.
- Positional GET stays byte-compatible with samtools on the verified path.
- `.fai` emission fails instead of describing non-uniform layout incorrectly.
- `.zfi` wire bytes remain versioned and backward-aware.
- Duplicate lookup is first exact match.
- Positional and BED coordinate systems remain distinct and explicit.
- Shared type detection changes are checked in stats, GET, and validate.
- GET and index hot paths stay bounded by fixed sequence payload buffers.
- Validate fix preserves biological sequence bytes and never overwrites input.

## Documentation update map

When CLI behavior changes, update:

- the relevant command page;
- [Command cheat sheet](Command-Cheat-Sheet);
- affected [Recipes](Recipes);
- [Troubleshooting](Troubleshooting) when diagnostics or recovery change;
- [FAQ](FAQ) when the user decision changes;
- [Limits and supported formats](Limits-and-Supported-Formats) when a boundary moves;
- `_Sidebar.md` when pages or command names change;
- the repository `README.md` for landing-page claims.

When format behavior changes, also update [Index formats](Index-Formats) and the internal `.zfi` wire contract.

The authored Wiki pages live in `wiki/` in the main repository. Edit and review that directory first. Publication to the GitHub Wiki is manual for now, so do not leave a correction only in the live Wiki.

> [!IMPORTANT]
> Examples are product claims. Verify command, status, stdout, stderr, and filesystem effects before publishing them.

## Writing rules

- Use short direct sentences and active voice.
- Keep prose unwrapped in Markdown source.
- Use fenced code blocks with language tags.
- Use tables only when comparison or repeated numeric structure earns one.
- Use no more than one or two GitHub alerts per page.
- Use alerts for user decisions, not decoration.
- Use `<details>` for secondary output or advanced explanations.
- Avoid emojis, decorative Unicode, and em-dashes.
- Link to one owner instead of duplicating deep format or implementation detail.

## Test documentation links

Before publishing Wiki changes:

1. Check every sidebar page exists.
2. Check every internal page link uses the published filename.
3. Preview alerts and collapsed sections on GitHub.
4. Run representative commands for changed examples.
5. Confirm the documented version.
6. Search for removed flags and stale limits.

## Report a code change honestly

State what succeeded, what failed or was redone, what was skipped, and whether the work is complete. Include the exact gates run.
