//! CLI entry and subcommand dispatcher for `index`, `get`, `stats`, and `validate`.
//!
//! Re-exports core modules for tests. Argument parsing, usage text, and routing live here;
//! indexing, extraction, stats, and validation logic live in the sibling modules.

const std = @import("std");
const builtin = @import("builtin");

// --- Module imports ---
pub const index_format = @import("index_format.zig");
pub const indexer = @import("indexer.zig");
pub const getter = @import("getter.zig");
pub const stats = @import("stats.zig");
pub const validator = @import("validator.zig");

// --- Re-exports for tests ---
pub const IndexRecord = index_format.IndexRecord;
pub const ZfiHeader = index_format.ZfiHeader;
pub const ZFI_MAGIC = index_format.ZFI_MAGIC;
pub const validateFasta = indexer.validateFasta;

const printErrorAndExit = index_format.printErrorAndExit;

const VERSION = "0.3.1";
const CHUNK_SIZE_FLAG = "--chunk-size";
const STRAND_AWARE_FLAG = "--strand-aware";
const STRAND_AWARE_ALIAS = "--honor-strand";
const RC_FLAG = "--rc";
const COMPLEMENT_ONLY_FLAG = "--complement-only";
const REVERSE_ONLY_FLAG = "--reverse-only";
const ANNOTATE_RC_FLAG = "--annotate-rc";

const USAGE =
    \\usage: z-fasta <command> [options]
    \\
    \\Commands:
    \\  index    Build a FASTA index (.zfi binary or .fai text)
    \\  get      Extract a sequence or region from a FASTA file
    \\  stats    Show statistics for a FASTA file
    \\  validate Check FASTA structure, alphabets, headers, and fix safe format issues
    \\
    \\General options:
    \\  --help       Show this help message
    \\  --version    Print version
    \\
    \\Index options:
    \\  --emit-fai   Output FAI to stdout when every record has fixed line geometry;
    \\               otherwise fails and directs you to default .zfi indexing
    \\  --no-dedup   Keep duplicate sequence names in the index (default: first wins
    \\               at index time). get resolves duplicate names to the last record.
    \\
    \\Get usage:
    \\  z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt]
    \\               [--strand-aware] [--summary] [--rc|--complement-only|--reverse-only]
    \\               [--annotate-rc] [--chunk-size N|-1] <region> [region ...]
    \\  Region formats: NAME, NAME:START-END, NAME:START-
    \\  --strand-aware  Respect BED column 6; '-' emits reverse-complement output
    \\                  (alias: --honor-strand)
    \\  --rc            Reverse-complement extracted sequence output
    \\  --complement-only  Complement extracted sequence output without reversing
    \\  --reverse-only  Reverse extracted sequence output without complementing
    \\  --annotate-rc   Annotate transformed headers (for example: reverse complement)
    \\  --chunk-size N  Process BED rows in batches of N (default: 4096). Use 1 only
    \\                  for debugging: one row per batch, high per-row arena overhead.
    \\  --chunk-size -1 Process all BED rows in one batch (512 MiB cap)
    \\  Positional regions: max 1024 per invocation. --names loads the whole file
    \\  up to 512 MiB (--chunk-size does not stream --names). BED streams by default.
    \\
    \\Stats options:
    \\  --index-only   Only show index-derived stats (no composition scan)
    \\
    \\Validate usage:
    \\  z-fasta validate [options] <file.fasta>
    \\  --strict                 Treat warnings as errors
    \\  --json                   Emit JSON Lines instead of text
    \\  --summary                With --json, emit one summary object
    \\  --fix -o <file.fasta>    Write a fixed FASTA to a separate output path
    \\  --fix-format-only        With --fix, proceed despite alphabet errors
    \\  --schema NAME            Header schema: uniprot or refseq
    \\  --custom-alphabet CHARS  Override nucleotide/protein alphabet checks
    \\  --max-header-len N       Warn on headers longer than N bytes (default: 1024)
    \\  Event list stops at 10000 issues (error without --fix). --fix still streams to -o.
    \\
    \\Examples:
    \\  z-fasta index genome.fa                  Create .zfi binary index
    \\  z-fasta index --emit-fai genome.fa       Output FAI to stdout
    \\  z-fasta get genome.fa chr1:1000-2000     Extract region
    \\  z-fasta get genome.fa chr1               Extract full sequence
    \\  z-fasta get genome.fa chr1:1000-2000 --rc
    \\  z-fasta get genome.fa --bed regions.bed  Extract BED regions
    \\  z-fasta get genome.fa --names ids.txt    Extract whole sequences from a file
    \\  z-fasta get genome.fa --bed regions.bed --strand-aware --summary
    \\  z-fasta stats genome.fa                  Full stats with composition
    \\  z-fasta stats --index-only genome.fa     Quick index-only stats
    \\  z-fasta validate genome.fa               Check FASTA validity
    \\  z-fasta validate --json --summary genome.fa
    \\
;

const ParseChunkSizeError = error{InvalidChunkSize};
const ParseTransformFlagsError = error{ConflictingTransformFlags};

const ParsedTransformFlags = struct {
    orientation: getter.Orientation,
    annotate_transform: bool,
};

fn parseChunkSizeValue(raw: []const u8) ParseChunkSizeError!usize {
    if (std.mem.eql(u8, raw, "-1")) return getter.chunk_size_all;

    const parsed = std.fmt.parseInt(usize, raw, 10) catch {
        return error.InvalidChunkSize;
    };
    if (parsed == 0) return error.InvalidChunkSize;
    return parsed;
}

fn chunkSizeEqualsValue(arg: []const u8) ?[]const u8 {
    const prefix = CHUNK_SIZE_FLAG ++ "=";
    if (!std.mem.startsWith(u8, arg, prefix)) return null;
    return arg[prefix.len..];
}

fn parseTransformFlags(rc: bool, complement_only: bool, reverse_only: bool, annotate_transform: bool) ParseTransformFlagsError!ParsedTransformFlags {
    var selected: usize = 0;
    if (rc) selected += 1;
    if (complement_only) selected += 1;
    if (reverse_only) selected += 1;
    if (selected > 1) return error.ConflictingTransformFlags;

    const orientation = if (rc)
        getter.Orientation.reverseComplement()
    else if (complement_only)
        getter.Orientation.complementOnly()
    else if (reverse_only)
        getter.Orientation.reverseOnly()
    else
        getter.Orientation{};

    return .{
        .orientation = orientation,
        .annotate_transform = annotate_transform,
    };
}

fn printUsageAndExit() noreturn {
    std.debug.print("{s}", .{USAGE});
    std.process.exit(1);
}

/// Positional argv tokens must not start with `-`. Known flags are matched earlier.
fn rejectUnknownOption(arg: []const u8) void {
    if (arg.len > 0 and arg[0] == '-') {
        printErrorAndExit("error: unknown option: {s}\n", .{arg});
    }
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch {};
    std.process.exit(0);
}

fn printVersionAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, "z-fasta " ++ VERSION ++ "\n") catch {};
    std.process.exit(0);
}

fn writeStdoutAll(bytes: []const u8) void {
    if (comptime builtin.os.tag == .linux) {
        var written: usize = 0;
        while (written < bytes.len) {
            const result = std.os.linux.write(1, bytes[written..].ptr, bytes.len - written);
            if (std.os.linux.errno(result) != .SUCCESS) return;
            if (result == 0) return;
            written += result;
        }
    } else {
        var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
        defer threaded.deinit();
        std.Io.File.writeStreamingAll(.stdout(), threaded.io(), bytes) catch {};
    }
}

fn printHelpFastAndExit() noreturn {
    writeStdoutAll(USAGE);
    std.process.exit(0);
}

fn printVersionFastAndExit() noreturn {
    writeStdoutAll("z-fasta " ++ VERSION ++ "\n");
    std.process.exit(0);
}

pub fn main(init: std.process.Init.Minimal) void {
    var args = std.process.Args.Iterator.initAllocator(init.args, std.heap.page_allocator) catch {
        std.debug.print("error: out of memory\n", .{});
        std.process.exit(1);
    };
    defer args.deinit();
    _ = args.skip();

    const cmd = args.next() orelse {
        printUsageAndExit();
    };

    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        printHelpFastAndExit();
    }
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) {
        printVersionFastAndExit();
    }

    var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    if (std.mem.eql(u8, cmd, "index")) {
        runIndex(io, init.environ, &args);
    } else if (std.mem.eql(u8, cmd, "get")) {
        runGetCmd(io, &args);
    } else if (std.mem.eql(u8, cmd, "stats")) {
        runStatsCmd(io, &args);
    } else if (std.mem.eql(u8, cmd, "validate")) {
        runValidateCmd(io, &args);
    } else {
        printUsageAndExit();
    }
}

// ============================================================================
// Subcommand: index
// ============================================================================

fn runIndex(io: std.Io, environ: std.process.Environ, args: *std.process.Args.Iterator) void {
    var emit_fai = false;
    var enable_dedup = true;
    var fasta_path: ?[]const u8 = null;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) {
            printVersionAndExit(io);
        } else if (std.mem.eql(u8, arg, "--emit-fai")) {
            emit_fai = true;
        } else if (std.mem.eql(u8, arg, "--no-dedup")) {
            enable_dedup = false;
        } else if (std.mem.eql(u8, arg, "--dedup")) {
            enable_dedup = true;
        } else {
            rejectUnknownOption(arg);
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printUsageAndExit();
    };

    indexer.runIndex(io, environ, path, .{
        .emit_fai = emit_fai,
        .enable_dedup = enable_dedup,
    }) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
        error.SourceOpenFailed => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
        error.SourceStatFailed => printErrorAndExit("error: failed to stat file: {s}\n", .{path}),
        error.EmptyFile => printErrorAndExit("error: file is empty: {s}\n", .{path}),
        error.NotFasta => printErrorAndExit("error: not a FASTA file: {s}\n", .{path}),
        error.HeaderTooLong => printErrorAndExit(
            "error: sequence name exceeds {d} bytes: {s}\n",
            .{ indexer.max_index_name_len, path },
        ),
        error.NonUniformFai => printErrorAndExit(
            "error: cannot emit .fai for non-uniform sequence layout; run 'z-fasta index' (default) to write .zfi\n",
            .{},
        ),
        error.NoValidSequences => printErrorAndExit(
            "error: no valid sequences found in: {s}\n",
            .{path},
        ),
        error.SourceReadFailed => printErrorAndExit("error: failed to read file: {s}\n", .{path}),
        error.NoUsableFaiSpool => printErrorAndExit(
            "error: no usable temporary directory for FAI spool\n",
            .{},
        ),
        error.FaiSpoolWriteFailed => printErrorAndExit(
            "error: failed to write temporary FAI spool\n",
            .{},
        ),
        error.FaiSpoolReadFailed => printErrorAndExit(
            "error: failed to read temporary FAI spool\n",
            .{},
        ),
        error.StdoutReplayFailed => printErrorAndExit("error: failed to replay FAI stdout\n", .{}),
        error.StdoutFlushFailed => printErrorAndExit("error: failed to flush FAI stdout\n", .{}),
        error.SourceChanged => printErrorAndExit(
            "error: source changed while indexing: {s}\n",
            .{path},
        ),
        error.OutputPathTooLong => printErrorAndExit("error: path too long\n", .{}),
        error.ZfiWriteFailed => printErrorAndExit("error: write failed\n", .{}),
        error.ZfiFinalizeFailed => printErrorAndExit(
            "error: failed to finalize index: {s}.zfi\n",
            .{path},
        ),
        error.ProcessingFailed => printErrorAndExit("error: processing failed\n", .{}),
        error.OutOfMemory => printErrorAndExit("error: processing failed\n", .{}),
        else => printErrorAndExit("error: processing failed\n", .{}),
    };
}

// ============================================================================
// Subcommand: get
// ============================================================================

fn runGetCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;
    var bed_path: ?[]const u8 = null;
    var names_path: ?[]const u8 = null;
    var honor_strand = false;
    var summary = false;
    var chunk_size: usize = 4_096;
    var rc = false;
    var complement_only = false;
    var reverse_only = false;
    var annotate_transform = false;
    // Static buffer: up to 1024 region strings without heap allocation.
    var region_buf: [1024][]const u8 = undefined;
    var region_count: usize = 0;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--bed")) {
            bed_path = args.next() orelse {
                printErrorAndExit("error: --bed requires a path or '-'\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--names")) {
            names_path = args.next() orelse {
                printErrorAndExit("error: --names requires a path\n", .{});
            };
        } else if (std.mem.eql(u8, arg, STRAND_AWARE_FLAG) or std.mem.eql(u8, arg, STRAND_AWARE_ALIAS)) {
            honor_strand = true;
        } else if (std.mem.eql(u8, arg, "--summary")) {
            summary = true;
        } else if (std.mem.eql(u8, arg, RC_FLAG)) {
            rc = true;
        } else if (std.mem.eql(u8, arg, COMPLEMENT_ONLY_FLAG)) {
            complement_only = true;
        } else if (std.mem.eql(u8, arg, REVERSE_ONLY_FLAG)) {
            reverse_only = true;
        } else if (std.mem.eql(u8, arg, ANNOTATE_RC_FLAG)) {
            annotate_transform = true;
        } else if (std.mem.eql(u8, arg, CHUNK_SIZE_FLAG)) {
            const raw = args.next() orelse {
                printErrorAndExit("error: --chunk-size requires a positive integer or -1\n", .{});
            };
            chunk_size = parseChunkSizeValue(raw) catch {
                printErrorAndExit("error: --chunk-size requires a positive integer or -1\n", .{});
            };
        } else if (chunkSizeEqualsValue(arg)) |raw| {
            chunk_size = parseChunkSizeValue(raw) catch {
                printErrorAndExit("error: --chunk-size requires a positive integer or -1\n", .{});
            };
        } else if (fasta_path == null) {
            rejectUnknownOption(arg);
            fasta_path = arg;
        } else {
            rejectUnknownOption(arg);
            if (region_count >= region_buf.len) {
                printErrorAndExit("error: too many regions (max 1024)\n", .{});
            }
            region_buf[region_count] = arg;
            region_count += 1;
        }
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt] [--strand-aware] [--summary] [--rc|--complement-only|--reverse-only] [--annotate-rc] [--chunk-size N|-1] <region> [region ...]\n", .{});
    };
    if (region_count == 0 and bed_path == null and names_path == null) {
        printErrorAndExit("error: usage: z-fasta get <file.fasta> [--bed file.bed|-] [--names file.txt] [--strand-aware] [--summary] [--rc|--complement-only|--reverse-only] [--annotate-rc] [--chunk-size N|-1] <region> [region ...]\n", .{});
    }

    const transform = parseTransformFlags(rc, complement_only, reverse_only, annotate_transform) catch {
        printErrorAndExit("error: --rc, --complement-only, and --reverse-only are mutually exclusive\n", .{});
    };
    if (transform.annotate_transform and transform.orientation.isIdentity()) {
        printErrorAndExit("error: --annotate-rc requires --rc, --complement-only, or --reverse-only\n", .{});
    }

    getter.runGetWithOptions(io, path, .{
        .region_strs = region_buf[0..region_count],
        .bed_path = bed_path,
        .names_path = names_path,
        .honor_strand = honor_strand,
        .summary = summary,
        .chunk_size = chunk_size,
        .orientation = transform.orientation,
        .annotate_transform = transform.annotate_transform,
    });
}

test "parseChunkSizeValue accepts positive sizes and all sentinel" {
    try std.testing.expectEqual(@as(usize, 4096), try parseChunkSizeValue("4096"));
    try std.testing.expectEqual(getter.chunk_size_all, try parseChunkSizeValue("-1"));
}

test "parseChunkSizeValue rejects zero and invalid input" {
    try std.testing.expectError(error.InvalidChunkSize, parseChunkSizeValue("0"));
    try std.testing.expectError(error.InvalidChunkSize, parseChunkSizeValue("abc"));
}

test "chunkSizeEqualsValue parses inline assignment syntax" {
    try std.testing.expectEqualStrings("-1", chunkSizeEqualsValue("--chunk-size=-1").?);
    try std.testing.expect(chunkSizeEqualsValue("--chunk-size") == null);
}

test "parseTransformFlags returns requested orientation" {
    const rc = try parseTransformFlags(true, false, false, false);
    try std.testing.expect(rc.orientation.reverse);
    try std.testing.expect(rc.orientation.complement);

    const complement_only = try parseTransformFlags(false, true, false, true);
    try std.testing.expect(!complement_only.orientation.reverse);
    try std.testing.expect(complement_only.orientation.complement);
    try std.testing.expect(complement_only.annotate_transform);

    const reverse_only = try parseTransformFlags(false, false, true, false);
    try std.testing.expect(reverse_only.orientation.reverse);
    try std.testing.expect(!reverse_only.orientation.complement);
}

test "parseTransformFlags rejects conflicting transform flags" {
    try std.testing.expectError(error.ConflictingTransformFlags, parseTransformFlags(true, true, false, false));
    try std.testing.expectError(error.ConflictingTransformFlags, parseTransformFlags(true, false, true, false));
    try std.testing.expectError(error.ConflictingTransformFlags, parseTransformFlags(false, true, true, false));
}

// ============================================================================
// Subcommand: stats
// ============================================================================

fn runStatsCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;
    var index_only = false;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--index-only")) {
            index_only = true;
        } else {
            rejectUnknownOption(arg);
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta stats [--index-only] <file.fasta>\n", .{});
    };

    stats.runStats(io, path, index_only);
}

// ============================================================================
// Subcommand: validate
// ============================================================================

fn runValidateCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;
    var options = validator.Options{};
    var saw_json = false;
    var saw_summary = false;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--strict")) {
            options.strict = true;
        } else if (std.mem.eql(u8, arg, "--json")) {
            saw_json = true;
            options.output_mode = if (saw_summary) .json_summary else .json_lines;
        } else if (std.mem.eql(u8, arg, "--summary")) {
            saw_summary = true;
            options.output_mode = .json_summary;
        } else if (std.mem.eql(u8, arg, "--fix")) {
            options.fix = true;
        } else if (std.mem.eql(u8, arg, "--fix-format-only")) {
            options.fix_format_only = true;
        } else if (std.mem.eql(u8, arg, "-o")) {
            options.output_path = args.next() orelse {
                printErrorAndExit("error: -o requires an output path\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--schema")) {
            const raw = args.next() orelse {
                printErrorAndExit("error: --schema requires uniprot or refseq\n", .{});
            };
            options.schema = parseValidateSchema(raw);
        } else if (std.mem.eql(u8, arg, "--custom-alphabet")) {
            options.custom_alphabet = args.next() orelse {
                printErrorAndExit("error: --custom-alphabet requires characters\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--max-header-len")) {
            const raw = args.next() orelse {
                printErrorAndExit("error: --max-header-len requires a positive integer\n", .{});
            };
            options.max_header_len = parsePositiveUsize(raw, "--max-header-len");
        } else if (fasta_path == null) {
            rejectUnknownOption(arg);
            fasta_path = arg;
        } else {
            rejectUnknownOption(arg);
            printErrorAndExit("error: validate accepts exactly one FASTA path\n", .{});
        }
    }

    if (saw_summary and !saw_json) {
        printErrorAndExit("error: validate --summary requires --json\n", .{});
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta validate [options] <file.fasta>\n", .{});
    };

    validator.runValidate(io, path, options);
}

fn parseValidateSchema(raw: []const u8) validator.Schema {
    if (std.mem.eql(u8, raw, "uniprot")) return .uniprot;
    if (std.mem.eql(u8, raw, "refseq")) return .refseq;
    printErrorAndExit("error: --schema must be uniprot or refseq\n", .{});
}

fn parsePositiveUsize(raw: []const u8, comptime flag: []const u8) usize {
    const parsed = std.fmt.parseInt(usize, raw, 10) catch {
        printErrorAndExit("error: {s} requires a positive integer\n", .{flag});
    };
    if (parsed == 0) {
        printErrorAndExit("error: {s} requires a positive integer\n", .{flag});
    }
    return parsed;
}
