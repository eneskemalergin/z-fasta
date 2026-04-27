const std = @import("std");
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

const VERSION = "0.2.3";

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
    \\  z-fasta get <file.fasta> <region> [region ...]
    \\  Region formats: NAME, NAME:START-END, NAME:START-
    \\
    \\Stats options:
    \\  --index-only   Only show index-derived stats (no composition scan)
    \\
    \\Examples:
    \\  z-fasta index genome.fa                  Create .zfi binary index
    \\  z-fasta index --emit-fai genome.fa       Output FAI to stdout
    \\  z-fasta get genome.fa chr1:1000-2000     Extract region
    \\  z-fasta get genome.fa chr1               Extract full sequence
    \\  z-fasta get genome.fa chr1 chr2 chrM     Extract multiple sequences
    \\  z-fasta stats genome.fa                  Full stats with composition
    \\  z-fasta stats --index-only genome.fa     Quick index-only stats
    \\
;

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

pub fn main(init: std.process.Init) void {
    const io = init.io;
    var args = std.process.Args.Iterator.init(init.minimal.args);
    _ = args.skip();

    const cmd = args.next() orelse {
        printUsageAndExit();
    };

    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        printHelpAndExit(io);
    }
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) {
        printVersionAndExit(io);
    }

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
    // Static buffer: up to 1024 region strings without heap allocation.
    var region_buf: [1024][]const u8 = undefined;
    var region_count: usize = 0;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
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
        printErrorAndExit("error: usage: z-fasta get <file.fasta> <region> [region ...]\n", .{});
    };
    if (region_count == 0) {
        printErrorAndExit("error: usage: z-fasta get <file.fasta> <region> [region ...]\n", .{});
    }

    getter.runGet(io, path, region_buf[0..region_count]);
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
