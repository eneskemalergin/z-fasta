const std = @import("std");
const builtin = @import("builtin");
const posix = std.posix;

// Module imports
pub const index_format = @import("index_format.zig");
pub const indexer = @import("indexer.zig");
pub const getter = @import("getter.zig");
pub const stats = @import("stats.zig");

// Re-exports for backward compatibility (tests import these)
pub const IndexRecord = index_format.IndexRecord;
pub const ZfiHeader = index_format.ZfiHeader;
pub const ZFI_MAGIC = index_format.ZFI_MAGIC;
pub const writeZfi = index_format.writeZfi;
pub const validateFasta = indexer.validateFasta;
pub const scanHeaders = indexer.scanHeaders;

const VERSION = "0.2.8";
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
    \\
    \\General options:
    \\  --help       Show this help message
    \\  --version    Print version
    \\
    \\Index options:
    \\  --emit-fai   Output FAI format to stdout (default: create .zfi file)
    \\  --no-dedup   Disable duplicate name filtering (default: dedup ON)
    \\  --low-mem    Use chunked reader instead of mmap (4 MB constant memory)
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
    \\  --chunk-size -1 Process all BED rows in one batch
    \\                  Default chunk size: 4096 BED rows
    \\
    \\Stats options:
    \\  --index-only   Only show index-derived stats (no composition scan)
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

fn printErrorAndExit(comptime fmt: []const u8, args_tuple: anytype) noreturn {
    std.debug.print(fmt, args_tuple);
    std.process.exit(1);
}

fn printUsageAndExit() noreturn {
    std.debug.print("{s}", .{USAGE});
    std.process.exit(1);
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
    var args = std.process.Args.Iterator.init(init.args);
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
        runIndex(io, &args);
    } else if (std.mem.eql(u8, cmd, "get")) {
        runGetCmd(io, &args);
    } else if (std.mem.eql(u8, cmd, "stats")) {
        runStatsCmd(io, &args);
    } else {
        printUsageAndExit();
    }
}

// ============================================================================
// Subcommand: index
// ============================================================================

fn runIndex(io: std.Io, args: *std.process.Args.Iterator) void {
    var emit_fai = false;
    var enable_dedup = true;
    var low_mem = false;
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
        } else if (std.mem.eql(u8, arg, "--low-mem")) {
            low_mem = true;
            emit_fai = true;
        } else {
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printUsageAndExit();
    };

    if (low_mem) {
        indexer.runChunkedMode(io, path);
        return;
    }

    // Standard mmap mode
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| {
        switch (err) {
            error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{path}),
            error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{path}),
            else => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
        }
    };
    defer file.close(io);

    const stat = file.stat(io) catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{path});
    };

    if (stat.size == 0) {
        printErrorAndExit("error: file is empty: {s}\n", .{path});
    }

    const data = posix.mmap(
        null,
        stat.size,
        .{ .READ = true },
        .{ .TYPE = .PRIVATE },
        file.handle,
        0,
    ) catch {
        printErrorAndExit("error: failed to mmap file: {s}\n", .{path});
    };
    defer posix.munmap(data);

    posix.madvise(data.ptr, data.len, posix.MADV.SEQUENTIAL) catch {};

    if (data.len == 0 or data[0] != '>') {
        printErrorAndExit("error: not a FASTA file: {s}\n", .{path});
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    if (emit_fai) {
        var out_buf: [65536]u8 = undefined;
        var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
        const record_count = indexer.streamingScan(data, &stdout_fw.interface, .fai, enable_dedup, arena.allocator()) catch {
            printErrorAndExit("error: failed to scan/write\n", .{});
        };
        stdout_fw.flush() catch {};

        if (record_count == 0) {
            printErrorAndExit("error: no valid sequences found in: {s}\n", .{path});
        }
    } else {
        var zfi_path_buf: [4096]u8 = undefined;
        const zfi_path = std.fmt.bufPrint(&zfi_path_buf, "{s}.zfi", .{path}) catch {
            printErrorAndExit("error: path too long\n", .{});
        };

        const out_file = std.Io.Dir.cwd().createFile(io, zfi_path, .{}) catch {
            printErrorAndExit("error: cannot create: {s}\n", .{zfi_path});
        };
        defer out_file.close(io);

        var file_buf: [65536]u8 = undefined;
        var file_fw = out_file.writer(io, &file_buf);
        const writer = &file_fw.interface;

        const dummy_header = index_format.ZfiHeader{
            .magic = index_format.ZFI_MAGIC,
            .record_count = 0,
            .source_size = data.len,
        };
        writer.writeAll(std.mem.asBytes(&dummy_header)) catch {
            printErrorAndExit("error: write failed\n", .{});
        };

        const record_count = indexer.streamingScan(data, writer, .zfi, enable_dedup, arena.allocator()) catch {
            printErrorAndExit("error: scan failed\n", .{});
        };

        file_fw.flush() catch {};

        if (record_count == 0) {
            std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
            printErrorAndExit("error: no valid sequences found in: {s}\n", .{path});
        }

        // Fix record_count in header
        file_fw.seekTo(4) catch {};
        file_fw.interface.writeInt(u32, record_count, .little) catch {};
        file_fw.flush() catch {};

        std.debug.print("wrote {s} ({d} sequences)\n", .{ zfi_path, record_count });
    }
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
            fasta_path = arg;
        } else {
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
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta stats [--index-only] <file.fasta>\n", .{});
    };

    stats.runStats(io, path, index_only);
}
