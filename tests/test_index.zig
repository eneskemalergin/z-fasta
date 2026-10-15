//! Index tests: scanning across read boundaries, sidecar loading and validation, gates fixtures.

const std = @import("std");
const builtin = @import("builtin");
const main = @import("main");

const support = @import("support.zig");

const diagnostics = main.index_format.diagnostics;

const ZFASTA_BIN = support.ZFASTA_BIN;
const expectCliFailure = support.expectCliFailure;
const expectCliResult = support.expectCliResult;
const expectUnknownOptionRejected = support.expectUnknownOptionRejected;
const uniqueArtifactPath = support.uniqueArtifactPath;
const readTestFile = support.readRequiredFile;
const writeFastaArtifact = support.writeFastaArtifact;
const statMtimeNs = support.statMtimeNs;
const scanFaiData = main.indexer.scanFaiData;
const INDEX_READ_BUFFER_BYTES = main.indexer.INDEX_READ_BUFFER_BYTES;
const IndexRecord = main.index_format.IndexRecord;
const ZfiHeader = main.index_format.ZfiHeader;
const ZFI_MAGIC = main.index_format.ZFI_MAGIC;
const c = std.c;
const io = std.testing.io;

fn zfiBytes(allocator: std.mem.Allocator, index: *const main.indexer.ZfiIndex, source_size: u64, source_mtime_ns: u64) ![]u8 {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    return main.indexer.zfiIndexToBytes(allocator, index, source_size, source_mtime_ns, &failure);
}

fn loadAndDeinitForAllocationCheck(allocator: std.mem.Allocator, fasta_path: []const u8, mode: main.index_format.LoadMode) !void {
    var idx = try support.loadIndex(allocator, io, fasta_path, mode);
    defer idx.deinit();
}

fn writeFastaAndRawZfi(allocator: std.mem.Allocator, stem: []const u8, fasta: []const u8, zfi_bytes: []const u8) !struct { fasta_path: []u8, zfi_path: []u8 } {
    const fasta_path = try writeFastaArtifact(allocator, stem, fasta);
    errdefer allocator.free(fasta_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    errdefer allocator.free(zfi_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const fasta_mtime = try statMtimeNs(fasta_path);

    var zfi_copy = try allocator.dupe(u8, zfi_bytes);
    defer allocator.free(zfi_copy);
    const footer_len = main.index_format.ZFI_NAME_FOOTER_BYTES;
    const id_len = main.index_format.ZFI_SOURCE_ID_BYTES;
    if (zfi_copy.len >= footer_len + id_len) {
        const id_off = zfi_copy.len - footer_len - id_len;
        if (std.mem.eql(u8, zfi_copy[id_off..][0..4], &main.index_format.ZFI_SOURCE_ID_MAGIC)) {
            std.mem.writeInt(u64, zfi_copy[id_off + 4 ..][0..8], fasta_mtime, .little);
        }
    }

    const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
    defer zfi_file.close(io);
    try std.Io.File.writeStreamingAll(zfi_file, io, zfi_copy);
    return .{ .fasta_path = fasta_path, .zfi_path = zfi_path };
}

fn writeFastaAndFai(allocator: std.mem.Allocator, stem: []const u8, fasta: []const u8, fai: []const u8) !struct { fasta_path: []u8, fai_path: []u8 } {
    const fasta_path = try writeFastaArtifact(allocator, stem, fasta);
    errdefer allocator.free(fasta_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    errdefer allocator.free(fai_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, fai);
    }
    return .{ .fasta_path = fasta_path, .fai_path = fai_path };
}

fn supportsPosixFutimens() bool {
    return switch (builtin.os.tag) {
        .linux, .dragonfly, .freebsd, .netbsd, .openbsd, .illumos, .macos, .ios, .tvos, .watchos, .visionos, .driverkit, .maccatalyst => true,
        else => false,
    };
}

fn uniqueArtifactDirPath(allocator: std.mem.Allocator, stem: []const u8) ![]u8 {
    const path = try uniqueArtifactPath(allocator, stem, "dir");
    errdefer allocator.free(path);
    try std.Io.Dir.cwd().createDirPath(io, path);
    return path;
}

fn expectDirectoryEmpty(path: []const u8) !void {
    var dir = try std.Io.Dir.cwd().openDir(io, path, .{ .iterate = true });
    defer dir.close(io);
    var entries = dir.iterate();
    try std.testing.expectEqual(null, try entries.next(io));
}

fn writeLargeIndexDescription(writer: *std.Io.Writer) !void {
    const padding = [_]u8{'x'} ** 4096;
    try writer.writeAll(">ignored ");
    for (0..32 * 1024 * 1024 / padding.len) |_| try writer.writeAll(&padding);
    try writer.writeAll("\n");
}

fn markFileStaleOneHourAgo(file: std.Io.File) !void {
    if (!supportsPosixFutimens()) return error.SkipZigTest;

    switch (builtin.os.tag) {
        .linux => {
            var now: std.os.linux.timespec = undefined;
            if (std.os.linux.errno(std.os.linux.clock_gettime(.REALTIME, &now)) != .SUCCESS) {
                return error.Unexpected;
            }
            const stale: std.os.linux.timespec = .{ .sec = now.sec - 3600, .nsec = now.nsec };
            if (std.os.linux.errno(std.os.linux.futimens(file.handle, &.{ stale, stale })) != .SUCCESS) {
                return error.Unexpected;
            }
        },
        .dragonfly,
        .freebsd,
        .netbsd,
        .openbsd,
        .illumos,
        .macos,
        .ios,
        .tvos,
        .watchos,
        .visionos,
        .driverkit,
        .maccatalyst,
        => {
            var now: c.timespec = undefined;
            if (c.clock_gettime(c.CLOCK.REALTIME, &now) != 0) return error.Unexpected;
            const stale: c.timespec = .{ .sec = now.sec - 3600, .nsec = now.nsec };
            if (c.futimens(file.handle, &.{ stale, stale }) != 0) return error.Unexpected;
        },
        else => return error.SkipZigTest,
    }
}

fn scanRecordsForTest(allocator: std.mem.Allocator, data: []const u8) !std.ArrayList(IndexRecord) {
    var records: std.ArrayList(IndexRecord) = .empty;
    errdefer records.deinit(allocator);

    const Sink = struct {
        records: *std.ArrayList(IndexRecord),
        allocator: std.mem.Allocator,

        fn emit(ctx: *@This(), emit_info: main.indexer.FastaRecordEmit) !void {
            try ctx.records.append(ctx.allocator, emit_info.record);
        }
    };

    var reader = std.Io.Reader.fixed(data);
    var read_buf: [4096]u8 = undefined;
    var sink = Sink{ .records = &records, .allocator = allocator };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    _ = try main.indexer.scanFastaReader(
        allocator,
        &reader,
        &read_buf,
        .{ .failure = &failure },
        &sink,
        Sink.emit,
    );
    return records;
}

fn scanZfiForAllocationCheck(allocator: std.mem.Allocator, data: []const u8) !void {
    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
}

fn scanFaiForAllocationCheck(allocator: std.mem.Allocator, data: []const u8) !void {
    var storage: [256]u8 = undefined;
    var output = std.Io.Writer.fixed(&storage);
    _ = try scanFaiData(allocator, data, &output);
}

fn fuzzFastaScanner(_: void, smith: *std.testing.Smith) !void {
    var storage: [512]u8 = undefined;
    const data = storage[0..smith.slice(&storage)];
    var index = try main.indexer.scanZfiData(std.testing.allocator, data);
    defer index.deinit(std.testing.allocator);

    for (index.records.items) |record| {
        const data_len: u64 = data.len;
        try std.testing.expect(record.seq_offset < data_len);
        try std.testing.expect(record.seq_len <= data_len - record.seq_offset);
        try std.testing.expect(record.line_bases > 0);
        try std.testing.expect(record.line_bytes >= record.line_bases);
        const name_offset = std.math.cast(usize, record.name_offset) orelse return error.TestUnexpectedResult;
        try std.testing.expect(name_offset <= index.name_blob.items.len);
        if (name_offset > index.name_blob.items.len) return error.TestUnexpectedResult;
        try std.testing.expect(record.name_len <= index.name_blob.items.len - name_offset);
    }
}

fn requireFixturePath(path: []const u8) !void {
    std.Io.Dir.cwd().access(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            std.debug.print("required test fixture missing: {s}\n", .{path});
            return error.RequiredFixtureMissing;
        },
        else => |e| return e,
    };
}

fn scanFaiReaderForTest(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    read_buf: []u8,
    writer: *std.Io.Writer,
) !u32 {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    return main.indexer.scanFaiReader(allocator, reader, read_buf, writer, .{ .failure = &failure });
}

fn expectReaderMatchesProduction(
    allocator: std.mem.Allocator,
    data: []const u8,
    read_buf: []u8,
    label: []const u8,
) !void {
    var baseline_index = try main.indexer.scanZfiData(allocator, data);
    defer baseline_index.deinit(allocator);

    var zfi_reader = std.Io.Reader.fixed(data);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var fragmented_index = try main.indexer.scanZfiReader(
        allocator,
        &zfi_reader,
        read_buf,
        .{ .failure = &failure },
    );
    defer fragmented_index.deinit(allocator);

    const baseline_bytes = try zfiBytes(allocator, &baseline_index, data.len, 0);
    defer allocator.free(baseline_bytes);
    const fragmented_bytes = try zfiBytes(allocator, &fragmented_index, data.len, 0);
    defer allocator.free(fragmented_bytes);

    std.testing.expectEqualSlices(u8, baseline_bytes, fragmented_bytes) catch {
        std.debug.print("ZFI reader-size mismatch: {s}, buffer={d}\n", .{ label, read_buf.len });
        return error.TestExpectedEqual;
    };

    var baseline_fai = std.Io.Writer.Allocating.init(allocator);
    defer baseline_fai.deinit();
    const baseline_count = scanFaiData(allocator, data, &baseline_fai.writer) catch |baseline_err| {
        var fragmented_fai = std.Io.Writer.Allocating.init(allocator);
        defer fragmented_fai.deinit();
        var fai_reader = std.Io.Reader.fixed(data);
        try std.testing.expectError(
            baseline_err,
            scanFaiReaderForTest(
                allocator,
                &fai_reader,
                read_buf,
                &fragmented_fai.writer,
            ),
        );
        return;
    };

    var fragmented_fai = std.Io.Writer.Allocating.init(allocator);
    defer fragmented_fai.deinit();
    var fai_reader = std.Io.Reader.fixed(data);
    const fragmented_count = try scanFaiReaderForTest(
        allocator,
        &fai_reader,
        read_buf,
        &fragmented_fai.writer,
    );

    try std.testing.expectEqual(baseline_count, fragmented_count);
    std.testing.expectEqualStrings(baseline_fai.written(), fragmented_fai.written()) catch {
        std.debug.print("FAI reader-size mismatch: {s}, buffer={d}\n", .{ label, read_buf.len });
        return error.TestExpectedEqual;
    };
}

fn appendZfiHeader(buf: *std.ArrayList(u8), allocator: std.mem.Allocator, record_count: u32, source_size: u64) !void {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = record_count,
        .source_size = source_size,
    };
    try buf.appendSlice(allocator, std.mem.asBytes(&header));
}

fn expectFaiLoaderModesAgree(allocator: std.mem.Allocator, fasta_path: []const u8) !void {
    var full = try support.loadIndex(allocator, io, fasta_path, .on_access);
    defer full.deinit();
    const names = try allocator.alloc([]const u8, full.records.len());
    defer allocator.free(names);
    for (names, 0..) |*name, i| name.* = full.recordName(i).?;
    var requested = try support.loadIndex(allocator, io, fasta_path, .{ .requested = names });
    defer requested.deinit();
    var stats_scan = try support.loadIndex(allocator, io, fasta_path, .stats_scan);
    defer stats_scan.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, full.source);
    try std.testing.expectEqual(full.source, requested.source);
    try std.testing.expectEqual(full.source, stats_scan.source);
    try std.testing.expectEqual(full.records.len(), requested.records.len());
    try std.testing.expectEqual(full.records.len(), stats_scan.records.len());
    try std.testing.expectEqual(@as(usize, 0), stats_scan.name_slices.len);

    for (0..full.records.len()) |i| {
        const rec = try recordAt(&full, i);
        const streamed = try recordAt(&requested, i);
        const scanned = try recordAt(&stats_scan, i);

        try std.testing.expectEqual(rec.seq_offset, streamed.seq_offset);
        try std.testing.expectEqual(rec.seq_len, streamed.seq_len);
        try std.testing.expectEqual(rec.line_bases, streamed.line_bases);
        try std.testing.expectEqual(rec.line_bytes, streamed.line_bytes);

        try std.testing.expectEqual(rec.seq_offset, scanned.seq_offset);
        try std.testing.expectEqual(rec.seq_len, scanned.seq_len);
        try std.testing.expectEqual(rec.line_bases, scanned.line_bases);
        try std.testing.expectEqual(rec.line_bytes, scanned.line_bytes);

        try std.testing.expectEqualStrings(full.recordName(i).?, requested.recordName(i).?);
        try std.testing.expectEqual(null, stats_scan.recordName(i));
    }
}

fn expectCorruptFaiFixture(allocator: std.mem.Allocator, stem: []const u8, code: diagnostics.Code) !void {
    const src_fa = try std.fmt.allocPrint(allocator, "{s}.fasta", .{stem});
    const src_fai = try std.fmt.allocPrint(allocator, "{s}.fasta.fai", .{stem});
    const fasta_data = try readTestFile(allocator, src_fa);
    defer allocator.free(fasta_data);
    const fai_data = try readTestFile(allocator, src_fai);
    defer allocator.free(fai_data);

    const paths = try writeFastaAndFai(allocator, "gates-fai", fasta_data, fai_data);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, code);
}

fn recordAt(index: *const main.index_format.LoadedIndex, i: usize) !main.index_format.IndexRecord {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    return index.checkedRecord(i, &failure);
}

fn expectLoadFailure(allocator: std.mem.Allocator, path: []const u8, mode: main.index_format.LoadMode, code: diagnostics.Code) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(if (mode == .stats_scan) .stats else .get);
    try std.testing.expectError(diagnostics.metadata(code).returned_error, main.index_format.loadIndex(allocator, io, path, .automatic, mode, &failure));
    try std.testing.expectEqual(code, failure.primary().?.code());
}

fn keepIndexNewer(index_path: []const u8, fasta_path: []const u8) !void {
    const later = std.Io.Timestamp.fromNanoseconds(@as(i96, @intCast(try statMtimeNs(fasta_path) + std.time.ns_per_s)));
    try std.Io.Dir.cwd().setTimestamps(io, index_path, .{ .modify_timestamp = .{ .new = later } });
}

fn expectIndexLoadCode(allocator: std.mem.Allocator, path: []const u8, mode: main.index_format.LoadMode, code: diagnostics.Code, failure: *diagnostics.Failure) !void {
    failure.reset(.get);
    if (main.index_format.loadIndex(allocator, io, path, .automatic, mode, failure)) |loaded| {
        var index = loaded;
        defer index.deinit();
        if (index.on_demand == null) return error.ExpectedIndexFailure;
        const err = for (0..index.records.len()) |i| {
            _ = index.checkedRecord(i, failure) catch |err| break err;
        } else return error.ExpectedIndexFailure;
        try std.testing.expectEqual(diagnostics.metadata(code).returned_error, err);
        try std.testing.expectEqual(code, failure.primary().?.code());
    } else |err| {
        try std.testing.expectEqual(diagnostics.metadata(code).returned_error, err);
        try std.testing.expectEqual(code, failure.primary().?.code());
        try std.testing.expect(!failure.complete);
    }
}

fn failIndexRead(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
    return error.InputOutput;
}

fn shortIndexRead(_: ?*anyopaque, _: std.Io.File, _: []const []u8, _: u64) std.Io.File.ReadPositionalError!usize {
    return 0;
}

fn failIndexStat(_: ?*anyopaque, _: std.Io.File) std.Io.File.StatError!std.Io.File.Stat {
    return error.AccessDenied;
}

fn failIndexCreate(_: ?*anyopaque, _: std.Io.Dir, _: []const u8, _: std.Io.Dir.CreateFileAtomicOptions) std.Io.Dir.CreateFileAtomicError!std.Io.File.Atomic {
    return error.AccessDenied;
}

fn failIndexWrite(_: ?*anyopaque, _: std.Io.File, _: []const u8, _: []const []const u8, _: usize, _: u64) std.Io.File.WritePositionalError!usize {
    return error.NoSpaceLeft;
}

fn failIndexPublish(_: ?*anyopaque, _: std.Io.Dir, _: []const u8, _: std.Io.Dir, _: []const u8) std.Io.Dir.RenameError!void {
    return error.AccessDenied;
}

fn disappearIndexSource(_: ?*anyopaque, _: std.Io.Dir, _: []const u8, _: std.Io.Dir.StatFileOptions) std.Io.Dir.StatFileError!std.Io.File.Stat {
    return error.FileNotFound;
}

fn failIndexSpoolCreate(_: ?*anyopaque, _: std.Io.Dir, _: []const u8, _: std.Io.Dir.CreateFileOptions) std.Io.File.OpenError!std.Io.File {
    return error.AccessDenied;
}

fn failIndexCleanup(_: ?*anyopaque, _: std.Io.Dir, _: []const u8) std.Io.Dir.DeleteFileError!void {
    return error.AccessDenied;
}

fn discardIndexStdout(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
    switch (operation) {
        .file_write_streaming => |write| {
            var bytes = write.header.len;
            for (write.data, 0..) |part, i| bytes += part.len * (if (i + 1 == write.data.len) write.splat else 1);
            return .{ .file_write_streaming = bytes };
        },
        else => return std.testing.io.vtable.operate(userdata, operation),
    }
}

fn runIndexForAllocationCheck(allocator: std.mem.Allocator, path: []const u8, emit_fai: bool) !void {
    var vtable = io.vtable.*;
    vtable.operate = discardIndexStdout;
    var failure: diagnostics.Failure = undefined;
    _ = main.indexer.runIndex(allocator, .{ .userdata = io.userdata, .vtable = &vtable }, .empty, path, .{ .emit_fai = emit_fai }, &failure) catch |err| {
        if (err == error.OutOfMemory) try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
        return err;
    };
}

fn failIndexStdout(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
    switch (operation) {
        .file_write_streaming => return .{ .file_write_streaming = error.BrokenPipe },
        else => return std.testing.io.vtable.operate(userdata, operation),
    }
}

fn expectIndexRunCode(faulty_io: std.Io, path: []const u8, emit_fai: bool, code: diagnostics.Code, failure: *diagnostics.Failure) !void {
    try std.testing.expectError(diagnostics.metadata(code).returned_error, main.indexer.runIndex(std.testing.allocator, faulty_io, .empty, path, .{ .emit_fai = emit_fai }, failure));
    try std.testing.expectEqual(code, failure.primary().?.code());
    try std.testing.expect(!failure.complete);
}

// --- FASTA reader tests ---
test "[cli] - [index]: rejects unknown options regardless of position" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";

    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", "--not-a-flag", fasta }, "--not-a-flag");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", fasta, "--not-a-flag" }, "--not-a-flag");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", "-z", fasta }, "-z");
    try expectUnknownOptionRejected(allocator, &.{ ZFASTA_BIN, "index", "--low-mem", fasta }, "--low-mem");
}

test "[cli] - [index]: rejects multiple and nonexistent FASTA paths" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "index", "tests/data/simple.fasta", "tests/data/single.fasta" },
        1,
        "error [usage.extra_argument] index: Too many arguments. argument=fasta_path allowed_count=1\nhint: Correct the command arguments; use --help for accepted options.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "index", "tests/data/definitely-missing-cli-failure.fasta" },
        1,
        "error [io.file_not_found] index: \"tests/data/definitely-missing-cli-failure.fasta\": Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n",
    );
    try expectCliFailure(
        allocator,
        &.{ ZFASTA_BIN, "index", "tests/data/simple.fasta/child" },
        1,
        "error [io.open_failed] index: \"tests/data/simple.fasta/child\": Could not open the required input.; cause=\"NotDir\"\nhint: Check the named input/output and retained I/O cause.\n",
    );
}

test "[cli] - [index]: rejects inputs without indexable records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const empty_records = try writeFastaArtifact(allocator, "index-empty-records", ">a\n>b\n");
    defer std.Io.Dir.cwd().deleteFile(io, empty_records) catch {};
    const Case = struct { path: []const u8, code: diagnostics.Code };
    for ([_]Case{
        .{ .path = "tests/data/empty.fasta", .code = .@"fasta.no_records" },
        .{ .path = "tests/data/not_fasta.txt", .code = .@"index.initial_header_required" },
        .{ .path = empty_records, .code = .@"index.no_indexable_records" },
    }) |case| {
        for ([_]bool{ false, true }) |fai| {
            const argv: []const []const u8 = if (fai) &.{ ZFASTA_BIN, "index", "--emit-fai", case.path } else &.{ ZFASTA_BIN, "index", case.path };
            try support.expectCliFailureCode(allocator, argv, case.code);
        }
    }
}

test "[cli] - [index]: atomically replaces ZFI and keeps duplicate names" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "index-cli-success", ">dup\nA\n>dup\nCC\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const tmp_path = try std.fmt.allocPrint(allocator, "{s}.tmp", .{zfi_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteTree(io, tmp_path) catch {};
    try std.Io.Dir.cwd().createDir(io, tmp_path, .default_dir);
    const all_message = try std.fmt.allocPrint(allocator, "wrote {s} (2 sequences)\n", .{zfi_path});
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "index", fasta_path }, 0, "", all_message);
    {
        var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
        defer idx.deinit();
        try std.testing.expectEqual(@as(usize, 2), idx.records.len());
        try std.testing.expectEqual(@as(u64, 1), (try recordAt(&idx, 0)).seq_len);
        try std.testing.expectEqual(@as(u64, 2), (try recordAt(&idx, 1)).seq_len);
    }
    try std.Io.Dir.cwd().access(io, tmp_path, .{});
    for ([_][]const u8{ "--dedup", "--no-dedup" }) |option| {
        try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "index", option, fasta_path }, .@"usage.unknown_option");
    }
}

test "[cli] - [index]: large catalogs preserve both formats and duplicate names" {
    const allocator = std.testing.allocator;
    for ([_]bool{ false, true }) |compressed| {
        for ([_]bool{ false, true }) |non_uniform| {
            const fasta_path = try uniqueArtifactPath(allocator, "index-large-catalog", "fa");
            defer allocator.free(fasta_path);
            defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
            const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
            defer allocator.free(zfi_path);
            defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
            var expected_all_fai = std.Io.Writer.Allocating.init(allocator);
            defer expected_all_fai.deinit();
            {
                const file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
                defer file.close(io);
                var buffer: [4096]u8 = undefined;
                var output = file.writer(io, &buffer);
                try writeLargeIndexDescription(&output.interface);
                for (0..2) |pass| {
                    for (0..128) |i| {
                        var name_buffer: [64]u8 = undefined;
                        const name = if (compressed)
                            try std.fmt.bufPrint(&name_buffer, "ENST000000{d:0>6}.1", .{i})
                        else
                            try std.fmt.bufPrint(&name_buffer, "name-{d}", .{i});
                        const seq_offset = output.logicalPos() + name.len + 2;
                        try output.interface.print(">{s}\n{s}", .{
                            name,
                            if (pass == 0) "ACGT\nAC\n" else "TTTT\nTT\n",
                        });
                        if (non_uniform and i % 2 != 0) try output.interface.writeAll("GGGG\n");
                        try expected_all_fai.writer.print("{s}\t6\t{d}\t4\t5\n", .{ name, seq_offset });
                    }
                }
                try output.flush();
            }
            const file = try std.Io.Dir.cwd().openFile(io, fasta_path, .{});
            defer file.close(io);
            const stat = try file.stat(io);
            var buffer: [4096]u8 = undefined;
            var reader = file.reader(io, &buffer);
            try reader.seekTo(0);
            var read_buffer: [4096]u8 = undefined;
            var failure: diagnostics.Failure = undefined;
            failure.reset(.index);
            var expected_index = try main.indexer.scanZfiReader(
                allocator,
                &reader.interface,
                &read_buffer,
                .{ .failure = &failure },
            );
            defer expected_index.deinit(allocator);
            const expected_zfi = try zfiBytes(
                allocator,
                &expected_index,
                stat.size,
                try main.index_format.timestampToNs(stat.mtime),
            );
            defer allocator.free(expected_zfi);
            const message = try std.fmt.allocPrint(allocator, "wrote {s} (256 sequences)\n", .{zfi_path});
            defer allocator.free(message);
            try expectCliResult(allocator, &.{ ZFASTA_BIN, "index", fasta_path }, 0, "", message);
            const actual_zfi = try readTestFile(allocator, zfi_path);
            defer allocator.free(actual_zfi);
            try std.testing.expectEqualSlices(u8, expected_zfi, actual_zfi);
            const fai_args: []const []const u8 = &.{ ZFASTA_BIN, "index", "--emit-fai", fasta_path };
            if (non_uniform) {
                try support.expectCliFailureCode(allocator, fai_args, .@"index.nonuniform_fai");
            } else {
                try expectCliResult(allocator, fai_args, 0, expected_all_fai.written(), "");
            }
        }
    }
}

test "[cli] - [index]: large late failures preserve the prior index and remove temporary files" {
    const allocator = std.testing.allocator;
    const source_dir_path = try uniqueArtifactDirPath(allocator, "index-large-failure");
    defer allocator.free(source_dir_path);
    defer std.Io.Dir.cwd().deleteTree(io, source_dir_path) catch {};
    const spool_dir_path = try uniqueArtifactDirPath(allocator, "index-large-spool");
    defer allocator.free(spool_dir_path);
    defer std.Io.Dir.cwd().deleteTree(io, spool_dir_path) catch {};
    const fasta_path = try std.fmt.allocPrint(allocator, "{s}/source.fa", .{source_dir_path});
    defer allocator.free(fasta_path);
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer allocator.free(zfi_path);
    {
        const file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer file.close(io);
        var buffer: [4096]u8 = undefined;
        var output = file.writer(io, &buffer);
        try writeLargeIndexDescription(&output.interface);
        for (0..2000) |i| try output.interface.print(">name-{d}\nACGT\n", .{i});
        try output.interface.writeByte('>');
        for (0..16) |_| try output.interface.writeAll(&([_]u8{'A'} ** 4096));
        try output.interface.writeAll("\nACGT\n");
        try output.flush();
    }
    {
        const previous = try std.Io.Dir.cwd().createFile(io, zfi_path, .{});
        defer previous.close(io);
        try previous.writePositionalAll(io, "previous index", 0);
    }
    var env = try std.testing.environ.createMap(allocator);
    defer env.deinit();
    try env.put("TMPDIR", spool_dir_path);
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    for ([_][]const []const u8{
        &.{ ZFASTA_BIN, "index", "--error-format=json", fasta_path },
        &.{ ZFASTA_BIN, "index", "--error-format=json", "--emit-fai", fasta_path },
    }) |argv| {
        const result = try std.process.run(allocator, threaded.io(), .{
            .argv = argv,
            .environ_map = &env,
            .stdout_limit = .limited(128 * 1024),
            .stderr_limit = .limited(64 * 1024),
        });
        defer allocator.free(result.stdout);
        defer allocator.free(result.stderr);
        try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, result.term);
        try std.testing.expectEqualStrings("", result.stdout);
        try support.expectDiagnostic(allocator, result.stderr, .@"index.identifier_too_long");
        const previous = try readTestFile(allocator, zfi_path);
        defer allocator.free(previous);
        try std.testing.expectEqualStrings("previous index", previous);
        try expectDirectoryEmpty(spool_dir_path);
        var source_dir = try std.Io.Dir.cwd().openDir(io, source_dir_path, .{ .iterate = true });
        defer source_dir.close(io);
        var entries = source_dir.iterate();
        var count: usize = 0;
        while (try entries.next(io)) |entry| {
            try std.testing.expect(std.mem.eql(u8, entry.name, "source.fa") or
                std.mem.eql(u8, entry.name, "source.fa.zfi"));
            count += 1;
        }
        try std.testing.expectEqual(@as(usize, 2), count);
    }
}

test "[unit] - [fasta scanner]: records names and fixed-width geometry" {
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanRecordsForTest(arena.allocator(), data);
    try std.testing.expectEqual(@as(usize, 2), records.items.len);

    try std.testing.expectEqualStrings("seq1", records.items[0].getName(data));
    try std.testing.expectEqualStrings("seq2", records.items[1].getName(data));
    try std.testing.expectEqual(@as(u64, 4), records.items[0].seq_len);
    try std.testing.expectEqual(@as(u64, 6), records.items[0].seq_offset);
    try std.testing.expectEqual(@as(u64, 17), records.items[1].seq_offset);
    try std.testing.expectEqual(@as(u32, 4), records.items[0].line_bases);
    try std.testing.expectEqual(@as(u32, 5), records.items[0].line_bytes);
}

test "[unit] - [fasta scanner]: retains identifiers and discards descriptions" {
    const data = ">myseq some description here\nACGT\n>sp|P12345|PROT_NAME OS=Human\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanRecordsForTest(arena.allocator(), data);
    try std.testing.expectEqual(@as(usize, 2), records.items.len);
    try std.testing.expectEqualStrings("myseq", records.items[0].getName(data));
    try std.testing.expectEqualStrings("sp|P12345|PROT_NAME", records.items[1].getName(data));
}

test "[property] - [fasta scanner]: preserves counts and geometry across line layouts" {
    const Case = struct {
        data: []const u8,
        seq_len: u64,
        line_bases: ?u32 = null,
        line_bytes: ?u32 = null,
    };
    const cases = [_]Case{
        .{ .data = ">seq1\nAAAA\nBBBB\nCCCC\n", .seq_len = 12, .line_bases = 4, .line_bytes = 5 },
        .{ .data = ">chrSynthetic\nAAAAAAAA\nCCCCCCCC\nGGGG\n", .seq_len = 20, .line_bases = 8, .line_bytes = 9 },
        .{ .data = ">seq1\nACGTACGT  \nACGT\n", .seq_len = 12 },
        .{ .data = ">a\nAAAA\n>next\n", .seq_len = 4 },
        .{ .data = ">seq\nACGT    \nACGTACGT\n", .seq_len = 12 },
        .{ .data = ">seq1\r\nACGT\r\n", .seq_len = 4, .line_bases = 4, .line_bytes = 6 },
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    for (cases) |case| {
        const records = try scanRecordsForTest(arena.allocator(), case.data);
        try std.testing.expectEqual(@as(usize, 1), records.items.len);
        try std.testing.expectEqual(case.seq_len, records.items[0].seq_len);
        if (case.line_bases) |expected| {
            try std.testing.expectEqual(expected, records.items[0].line_bases);
        }
        if (case.line_bytes) |expected| {
            try std.testing.expectEqual(expected, records.items[0].line_bytes);
        }
    }
}

test "[edge] - [fasta scanner]: skips empty sequences" {
    const data = ">empty\n>next\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanRecordsForTest(arena.allocator(), data);
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqualStrings("next", records.items[0].getName(data));
}

test "[edge] - [fasta scanner]: keeps every duplicate name" {
    const data = ">dup\nAAAA\n>other\nCCCC\n>dup\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const records = try scanRecordsForTest(arena.allocator(), data);
    try std.testing.expectEqual(@as(usize, 3), records.items.len);
    try std.testing.expectEqualStrings("dup", records.items[0].getName(data));
    try std.testing.expectEqualStrings("other", records.items[1].getName(data));
    try std.testing.expectEqualStrings("dup", records.items[2].getName(data));
    try std.testing.expectEqual(@as(u64, 27), records.items[2].seq_offset);
}

test "[edge] - [fai scanner]: keeps duplicate and distinct names in order" {
    const data = ">dup\nAAAA\n>other\nCCCC\n>dup\nGGGG\n>also\nTTTT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    const count = try scanFaiData(allocator, data, &fai.writer);

    try std.testing.expectEqual(@as(u32, 4), count);
    try std.testing.expectEqualStrings(
        "dup\t4\t5\t4\t5\nother\t4\t17\t4\t5\ndup\t4\t27\t4\t5\nalso\t4\t38\t4\t5\n",
        fai.written(),
    );
}

test "[failure] - [fasta scanner]: every allocation failure releases partial scans in both formats" {
    const data = ">alpha\nAAAA\n>beta\nCCCC\n>alpha\nGGGG\n";

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        scanZfiForAllocationCheck,
        .{data},
    );
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        scanFaiForAllocationCheck,
        .{data},
    );
}

test "[fuzz] - [fasta scanner]: returns bounded records for arbitrary bytes" {
    try std.testing.fuzz({}, fuzzFastaScanner, .{ .corpus = &.{
        "",
        ">a\nA\n",
        ">a description\r\nAC\r\nGT",
        ">a\nAAA\nCC\nGGGG\n",
        ">bad\xffname\nA\x00C\n",
    } });
}

// --- ZFI wire format ---
test "[integration] - [zfi writer]: creates the production layout" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);

    const path = try uniqueArtifactPath(allocator, "test-write", "zfi");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    try main.indexer.writeZfiIndexFile(io, path, &index, fasta_data.len, 0);

    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    const stat = try file.stat(io);

    var header_bytes: [@sizeOf(ZfiHeader)]u8 = undefined;
    const header_bytes_read = try std.Io.File.readPositionalAll(file, io, &header_bytes, 0);
    try std.testing.expectEqual(header_bytes.len, header_bytes_read);
    const header: ZfiHeader = @bitCast(header_bytes);

    try std.testing.expectEqualSlices(u8, &ZFI_MAGIC, &header.magic);
    try std.testing.expectEqual(@as(u32, 1), header.record_count);
    try std.testing.expectEqual(@as(u64, fasta_data.len), header.source_size);
    try std.testing.expectEqualStrings("seq", index.name_blob.items);
    try std.testing.expectEqual(
        @as(u64, 16 + 24 + index.name_blob.items.len + 4 * main.index_format.NameTable.slotCount(1) + main.index_format.ZFI_NAME_FOOTER_BYTES + main.index_format.ZFI_SOURCE_ID_BYTES),
        stat.size,
    );
}

test "[unit] - [zfi format]: retains frozen wire sizes" {
    try std.testing.expectEqual(@as(usize, 16), @sizeOf(ZfiHeader));
    try std.testing.expectEqual(@as(usize, 24), @sizeOf(main.index_format.ZfiRecord));
    try std.testing.expectEqual(@as(usize, 24), @sizeOf(main.index_format.SideTableLine));
}

test "[failure] - [zfi format]: rejects a pre-epoch timestamp" {
    try std.testing.expectEqual(
        @as(u64, 0),
        try main.index_format.timestampToNs(std.Io.Timestamp.fromNanoseconds(0)),
    );
    try std.testing.expectError(
        error.TimestampOutOfRange,
        main.index_format.timestampToNs(std.Io.Timestamp.fromNanoseconds(-1)),
    );
}

test "[property] - [zfi format]: encoders match the frozen little-endian layout" {
    const header = ZfiHeader{
        .magic = ZFI_MAGIC,
        .record_count = 0x01020304,
        .source_size = 0x0807060504030201,
    };
    try std.testing.expectEqualSlices(u8, &(ZFI_MAGIC ++ [_]u8{ 0x04, 0x03, 0x02, 0x01, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 }), &main.index_format.encodeZfiHeader(header));

    const ZfiRecord = main.index_format.ZfiRecord;
    var rec = IndexRecord{
        .name_offset = 0x1112131415,
        .name_len = 0x2122,
        .seq_offset = 0x313233343536,
        .seq_len = 0x414243444546,
        .line_bases = 0x51525354,
        .line_bytes = 0x51525356,
    };
    const uniform = ZfiRecord.encode(rec);
    try std.testing.expectEqual([3]u64{ 0x2122_313233343536, 0x1415_414243444546, 0x02_51525354_111213 }, uniform.words);
    try std.testing.expectEqual(rec, uniform.decode());

    try rec.markNonUniform(main.index_format.MAX_SIDE_TABLE_OFFSET);
    const side = ZfiRecord.encode(rec);
    try std.testing.expectEqual(@as(u64, 0xFFFFFFFFFF_111213), side.words[2]);
    const decoded = side.decode();
    try std.testing.expect(!decoded.isUniformWidth());
    try std.testing.expectEqual(@as(u64, main.index_format.MAX_SIDE_TABLE_OFFSET), decoded.sideTableOffset());
    try std.testing.expectEqual(rec.name_offset, decoded.name_offset);
    try std.testing.expectEqual(rec.seq_len, decoded.seq_len);

    const footer = main.index_format.encodeZfiNameFooter(0x0102030405060708);
    try std.testing.expectEqualSlices(u8, &main.index_format.ZFI_NAME_FOOTER_MAGIC, footer[0..4]);
    try std.testing.expectEqual(@as(u64, 0x0102030405060708), std.mem.readInt(u64, footer[4..12], .little));
}

test "[property] - [zfi writer]: emits encoded records and side-table counts" {
    const data = ">seq\nAAA\nCCCC\nGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_index = try main.indexer.scanZfiData(allocator, data);
    defer zfi_index.deinit(allocator);
    const bytes = try zfiBytes(allocator, &zfi_index, data.len, 0);
    defer allocator.free(bytes);

    const expected_header = main.index_format.encodeZfiHeader(.{
        .magic = ZFI_MAGIC,
        .record_count = 1,
        .source_size = data.len,
    });
    try std.testing.expectEqualSlices(u8, &expected_header, bytes[0..16]);

    const rec = zfi_index.records.items[0];
    try std.testing.expectEqualSlices(
        u8,
        std.mem.asBytes(&main.index_format.ZfiRecord.encode(rec)),
        bytes[16 .. 16 + 24],
    );

    const side_off: usize = @intCast(rec.sideTableOffset());
    try std.testing.expectEqual(
        @as(u64, 3),
        std.mem.readInt(u64, bytes[side_off..][0..8], .little),
    );
}

test "[unit] - [zfi scanner]: stores side tables for non-uniform records" {
    const data = ">seq\nAAA\nCCCC\nGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_index = try main.indexer.scanZfiData(allocator, data);
    defer zfi_index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 1), zfi_index.records.items.len);
    const rec = zfi_index.records.items[0];
    try std.testing.expect(!rec.isUniformWidth());
    try std.testing.expectEqual(
        @as(u64, @sizeOf(ZfiHeader) + main.index_format.ZFI_INDEX_RECORD_BYTES),
        rec.sideTableOffset(),
    );

    try std.testing.expectEqual(
        @as(u64, 3),
        std.mem.readInt(u64, zfi_index.side_tables.items[0..@sizeOf(u64)], .little),
    );
}

test "[integration] - [zfi loader]: resolves embedded record names" {
    const data = @embedFile("data/simple.fasta");
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    const zfi_bytes = try zfiBytes(allocator, &index, data.len, 0);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "stream-zfi", data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .on_access);
    defer idx.deinit();
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("seq1"));
    try std.testing.expectEqual(@as(?usize, 1), idx.lookupName("seq2"));
}

test "[edge] - [fasta scanner]: a sequence line longer than the read buffer keeps its geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', INDEX_READ_BUFFER_BYTES + 10);
    try fasta.appendSlice(allocator, "\n>seq2\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    const count = try scanFaiData(allocator, fasta.items, &output.writer);
    const expected_fai = try std.fmt.allocPrint(
        allocator,
        "seq1\t{d}\t6\t{d}\t{d}\nseq2\t4\t{d}\t4\t5\n",
        .{
            INDEX_READ_BUFFER_BYTES + 10,
            INDEX_READ_BUFFER_BYTES + 10,
            INDEX_READ_BUFFER_BYTES + 11,
            INDEX_READ_BUFFER_BYTES + 23,
        },
    );

    try std.testing.expectEqual(@as(u32, 2), count);
    try std.testing.expectEqualStrings(expected_fai, output.written());

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expectEqual(@as(u64, INDEX_READ_BUFFER_BYTES + 10), index.records.items[0].seq_len);
    try std.testing.expectEqual(@as(u64, 6), index.records.items[0].seq_offset);
    try std.testing.expectEqual(@as(u32, INDEX_READ_BUFFER_BYTES + 10), index.records.items[0].line_bases);
    try std.testing.expectEqual(@as(u32, INDEX_READ_BUFFER_BYTES + 11), index.records.items[0].line_bytes);
    try std.testing.expectEqual(@as(u64, 4), index.records.items[1].seq_len);
    try std.testing.expectEqual(@as(u64, INDEX_READ_BUFFER_BYTES + 23), index.records.items[1].seq_offset);
    const second_name = index.name_blob.items[index.records.items[1].name_offset..][0..index.records.items[1].name_len];
    try std.testing.expectEqualStrings("seq2", second_name);
}

test "[edge] - [fasta scanner]: preserves a header split at the read boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const chunk = INDEX_READ_BUFFER_BYTES;
    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', chunk - 8);
    try fasta.appendSlice(allocator, "\n>seq2\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    const count = try scanFaiData(allocator, fasta.items, &output.writer);
    const expected_fai = try std.fmt.allocPrint(
        allocator,
        "seq1\t{d}\t6\t{d}\t{d}\nseq2\t4\t{d}\t4\t5\n",
        .{ chunk - 8, chunk - 8, chunk - 7, chunk + 5 },
    );

    try std.testing.expectEqual(@as(u32, 2), count);
    try std.testing.expectEqualStrings(expected_fai, output.written());

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expectEqual(@as(u64, chunk - 8), index.records.items[0].seq_len);
    try std.testing.expectEqual(@as(u64, chunk + 5), index.records.items[1].seq_offset);
}

test "[property] - [fasta scanner]: reader sizes preserve both formats across token boundaries" {
    const data = ">alpha\tlong description\r\nACGT\r\nAC\r\n>beta other text\nTTAA";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const vector_len = std.simd.suggestVectorLength(u8) orelse 1;
    const line_bytes = 6;
    const read_sizes = [_]usize{
        1,
        2,
        if (vector_len > 1) vector_len - 1 else 1,
        vector_len,
        vector_len + 1,
        line_bytes - 1,
        line_bytes,
        line_bytes + 1,
        INDEX_READ_BUFFER_BYTES,
    };
    const read_storage = try allocator.alloc(u8, INDEX_READ_BUFFER_BYTES);
    {
        for (read_sizes) |read_size| {
            try expectReaderMatchesProduction(
                allocator,
                data,
                read_storage[0..read_size],
                "token-boundary-read-sizes",
            );
        }
    }
}

test "[property] - [fasta scanner]: enforces the header limit across read fragments" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', main.indexer.MAX_INDEX_NAME_LEN);
    try fasta.appendSlice(allocator, "\nACGT\n");

    var read_storage: [31]u8 = undefined;
    {
        try expectReaderMatchesProduction(
            allocator,
            fasta.items,
            &read_storage,
            "maximum-fragmented-header",
        );
    }

    try fasta.insert(allocator, main.indexer.MAX_INDEX_NAME_LEN + 1, 'A');
    {
        var zfi_reader = std.Io.Reader.fixed(fasta.items);
        var failure: diagnostics.Failure = undefined;
        failure.reset(.index);
        try std.testing.expectError(error.InvalidIndex, main.indexer.scanZfiReader(allocator, &zfi_reader, &read_storage, .{ .failure = &failure }));
        try std.testing.expectEqual(diagnostics.Code.@"index.identifier_too_long", failure.primary().?.code());

        var fai = std.Io.Writer.Allocating.init(allocator);
        defer fai.deinit();
        var fai_reader = std.Io.Reader.fixed(fasta.items);
        try std.testing.expectError(
            error.InvalidIndex,
            scanFaiReaderForTest(
                allocator,
                &fai_reader,
                &read_storage,
                &fai.writer,
            ),
        );
    }
}

test "[edge] - [fasta scanner]: preserves sequence data split at the read boundary" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq1\n");
    try fasta.appendNTimes(allocator, 'A', INDEX_READ_BUFFER_BYTES);
    try fasta.appendSlice(allocator, "ACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    const count = try scanFaiData(allocator, fasta.items, &output.writer);
    const expected = try std.fmt.allocPrint(
        allocator,
        "seq1\t{d}\t6\t{d}\t{d}\n",
        .{ INDEX_READ_BUFFER_BYTES + 4, INDEX_READ_BUFFER_BYTES + 4, INDEX_READ_BUFFER_BYTES + 5 },
    );

    try std.testing.expectEqual(@as(u32, 1), count);
    try std.testing.expectEqualStrings(expected, output.written());

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expectEqual(@as(u64, INDEX_READ_BUFFER_BYTES + 4), index.records.items[0].seq_len);
    try std.testing.expectEqual(@as(u32, INDEX_READ_BUFFER_BYTES + 4), index.records.items[0].line_bases);
    try std.testing.expectEqual(@as(u32, INDEX_READ_BUFFER_BYTES + 5), index.records.items[0].line_bytes);
}

test "[edge] - [fai scanner]: counts a wrapped short final line with newline" {
    const data = ">chrSynthetic\nAAAAAAAA\nCCCCCCCC\nGGGG\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var out = std.Io.Writer.Allocating.init(allocator);
    defer out.deinit();
    const record_count = try scanFaiData(allocator, data, &out.writer);

    try std.testing.expectEqual(@as(u32, 1), record_count);
    try std.testing.expectEqualStrings("chrSynthetic\t20\t14\t8\t9\n", out.written());
}

test "[failure] - [fasta scanner]: rejects names above the u16 limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', main.indexer.MAX_INDEX_NAME_LEN + 1);
    try fasta.appendSlice(allocator, "\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();

    try std.testing.expectError(error.InvalidIndex, scanFaiData(allocator, fasta.items, &output.writer));
}

test "[edge] - [fasta scanner]: starts the sequence at the first base-bearing line after a blank" {
    const data = ">seq\n\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    _ = try scanFaiData(allocator, data, &output.writer);
    try std.testing.expectEqualStrings("seq\t4\t6\t4\t5\n", output.written());

    var records = try scanRecordsForTest(allocator, data);
    defer records.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), records.items.len);
    try std.testing.expectEqual(@as(u64, 6), records.items[0].seq_offset);
    try std.testing.expectEqual(@as(u32, 4), records.items[0].line_bases);
    try std.testing.expectEqual(@as(u32, 5), records.items[0].line_bytes);
}

test "[edge] - [fai scanner]: accepts a short final line without terminal newline" {
    const data = ">seq\nAAAA\nBB";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    _ = try scanFaiData(allocator, data, &output.writer);
    try std.testing.expectEqualStrings("seq\t6\t5\t4\t5\n", output.written());
}

test "[failure] - [fasta scanner]: rejects trailing space as non-uniform for FAI" {
    const data = ">s\nAAAA \n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(!index.records.items[0].isUniformWidth());

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    try std.testing.expectError(
        error.InvalidIndex,
        scanFaiData(allocator, data, &fai.writer),
    );
}

test "[failure] - [fai scanner]: rejects a trailing tab as non-uniform" {
    const data = ">s\nAAAA\t\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var stream_fai = std.Io.Writer.Allocating.init(allocator);
    defer stream_fai.deinit();
    try std.testing.expectError(
        error.InvalidIndex,
        scanFaiData(allocator, data, &stream_fai.writer),
    );
}

test "[integration] - [fasta scanner]: preserves whitespace before the next header" {
    const allocator = std.testing.allocator;
    const cases = [_][]const u8{
        ">a\nAAAA\nAAAA\nAA A\n>b x\nTTTT\n",
        ">a\r\nAAAA\r\nAAAA\r\nAA A\r\n>b x\r\nTTTT\r\n",
        ">a\nAAAA\nAAAA\nAAA \n>b x\nTTTT\n",
    };
    for (cases) |data| {
        var index = try main.indexer.scanZfiData(allocator, data);
        defer index.deinit(allocator);
        try std.testing.expectEqual(@as(u64, 11), index.records.items[0].seq_len);
        try std.testing.expect(!index.records.items[0].isUniformWidth());
        const bytes = try zfiBytes(allocator, &index, data.len, 0);
        defer allocator.free(bytes);
        const paths = try writeFastaAndRawZfi(allocator, "stride-whitespace", data, bytes);
        defer allocator.free(paths.fasta_path);
        defer allocator.free(paths.zfi_path);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        var loaded = try support.loadIndex(allocator, io, paths.fasta_path, .on_access);
        defer loaded.deinit();
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        var output: [128]u8 = undefined;
        var writer = std.Io.Writer.fixed(&output);
        try main.getter.extractRegion(&loaded, "a:9-11", &writer, &failure);
        try std.testing.expectEqualStrings(">a:9-11\nAAA\n", writer.buffered());
        var fai = std.Io.Writer.fixed(&output);
        try std.testing.expectError(error.InvalidIndex, scanFaiData(allocator, data, &fai));
        try std.testing.expectEqual(@as(usize, 0), fai.buffered().len);
    }
}

test "[edge] - [fasta scanner]: accepts a final CRLF after a uniform LF body" {
    const data = ">seq\nAAAA\nAAAA\r\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    _ = try scanFaiData(allocator, data, &fai.writer);
    try std.testing.expectEqualStrings("seq\t8\t5\t4\t5\n", fai.written());
}

test "[failure] - [fasta scanner]: rejects mixed interior LF and CRLF as non-uniform for FAI" {
    const data = ">seq\nAAAA\nAAAA\r\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    try std.testing.expectError(
        error.InvalidIndex,
        scanFaiData(allocator, data, &fai.writer),
    );
}

test "[edge] - [zfi scanner]: marks interior blank lines non-uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const data = ">seq\nAAAA\n\nBBBB\n";

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
}

test "[edge] - [fai scanner]: skips empty records" {
    const data = ">empty\n>next\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    _ = try scanFaiData(allocator, data, &output.writer);
    try std.testing.expectEqualStrings("next\t4\t13\t4\t5\n", output.written());
}

test "[integration] - [index loader]: resolves an empty sequence name in both formats" {
    const data = ">\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    _ = try scanFaiData(allocator, data, &fai.writer);
    try std.testing.expectEqualStrings("\t4\t2\t4\t5\n", fai.written());

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expectEqual(@as(u16, 0), index.records.items[0].name_len);

    const zfi_bytes = try zfiBytes(allocator, &index, data.len, 0);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "empty-name-zfi", data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .on_access);
    defer idx.deinit();
    try std.testing.expectEqual(@as(usize, 1), idx.records.len());
    try std.testing.expectEqualStrings("", idx.recordName(0).?);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName(""));
}

test "[edge] - [fasta scanner]: keeps every empty sequence name" {
    const data = ">\nACGT\n>\nTTAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    _ = try scanFaiData(allocator, data, &fai.writer);
    try std.testing.expectEqualStrings("\t4\t2\t4\t5\n\t4\t9\t4\t5\n", fai.written());

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expectEqual(@as(u64, 4), index.records.items[1].seq_len);
}

test "[integration] - [fai loader]: resolves an empty sequence name" {
    const data = ">\nACGT\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "empty-name-fai", data, "\t4\t2\t4\t5\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer idx.deinit();
    try std.testing.expectEqual(@as(usize, 1), idx.records.len());
    try std.testing.expectEqualStrings("", idx.recordName(0).?);
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName(""));
}

test "[edge] - [fai scanner]: accepts long names within the limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.append(allocator, '>');
    try fasta.appendNTimes(allocator, 'A', 10_000);
    try fasta.appendSlice(allocator, " description\nACGT\n");

    var output = std.Io.Writer.Allocating.init(allocator);
    defer output.deinit();
    const count = try scanFaiData(allocator, fasta.items, &output.writer);
    try std.testing.expectEqual(@as(u32, 1), count);
    const suffix = "\t4\t10014\t4\t5\n";
    try std.testing.expectEqual(@as(usize, 10_000 + suffix.len), output.written().len);
    try std.testing.expectEqualSlices(u8, fasta.items[1..10_001], output.written()[0..10_000]);
    try std.testing.expectEqualStrings(suffix, output.written()[10_000..]);
}

test "[failure] - [zfi loader]: rejects record names outside the name blob" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try writeFastaArtifact(allocator, "corrupt-zfi", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    index.records.items[0].name_offset = 999_999;
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, 11, try statMtimeNs(fasta_path));

    try expectLoadFailure(std.testing.allocator, fasta_path, .catalog, .@"index.name_out_of_bounds");
}

test "[failure] - [zfi loader]: rejects a sequence offset at EOF" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try writeFastaArtifact(allocator, "past-end-zfi", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    index.records.items[0].seq_offset = fasta_data.len;
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    try expectLoadFailure(std.testing.allocator, fasta_path, .catalog, .@"index.sequence_offset_out_of_bounds");
}

test "[failure] - [zfi loader]: rejects a zero record count" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayList(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 0, 11);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-zero-count", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.no_records");
}

test "[failure] - [zfi loader]: rejects a truncated header" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndRawZfi(allocator, "zfi-trunc-header", ">seq1\nACGT\n", ZFI_MAGIC[0..3]);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.truncated");
}

test "[failure] - [zfi loader]: rejects a record count beyond the file" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayList(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 2, 11);
    const rec = IndexRecord{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&rec));

    const paths = try writeFastaAndRawZfi(allocator, "zfi-count-overflow", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.truncated");
}

test "[failure] - [zfi loader]: rejects invalid magic" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayList(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    var header = ZfiHeader{ .magic = .{ 'X', 'F', 'I', 0x01 }, .record_count = 1, .source_size = 11 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&header));
    const rec = IndexRecord{ .name_offset = 1, .name_len = 4, .seq_offset = 6, .seq_len = 4, .line_bases = 4, .line_bytes = 5 };
    try zfi_buf.appendSlice(allocator, std.mem.asBytes(&rec));

    const paths = try writeFastaAndRawZfi(allocator, "zfi-bad-magic", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.invalid_magic");
}

test "[failure] - [zfi loader]: rejects a name blob overlapping records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    var zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    defer allocator.free(zfi_bytes);

    const footer = main.index_format.encodeZfiNameFooter(zfi_bytes.len);
    @memcpy(zfi_bytes[zfi_bytes.len - footer.len ..], &footer);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-blob-overlap", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.name_blob_out_of_bounds");
}

test "[failure] - [zfi loader]: rejects a side table offset inside records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());

    try index.records.items[0].markNonUniform(@sizeOf(ZfiHeader));

    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-side-in-records", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, .@"index.side_table_offset_out_of_bounds");
}

test "[failure] - [zfi loader]: rejects an empty side table" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);

    const offset: usize = @intCast(index.records.items[0].sideTableOffset());
    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    defer allocator.free(zfi_bytes);
    @memset(zfi_bytes[offset .. offset + @sizeOf(u64)], 0);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-side-zero-lines", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, .@"index.side_table_empty");
}

test "[failure] - [zfi loader]: rejects malformed side-table lines" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    const table_offset: usize = @intCast(index.records.items[0].sideTableOffset());

    const Mutation = enum {
        broken_base_prefix,
        repeated_byte_offset,
        overlapping_line_span,
        line_past_fasta,
        wrong_base_total,
    };
    const cases = [_]struct { Mutation, diagnostics.Code }{
        .{ .broken_base_prefix, .@"index.side_base_start_mismatch" },
        .{ .repeated_byte_offset, .@"index.side_line_order_invalid" },
        .{ .overlapping_line_span, .@"index.side_line_order_invalid" },
        .{ .line_past_fasta, .@"index.sequence_span_out_of_bounds" },
        .{ .wrong_base_total, .@"index.side_sequence_length_mismatch" },
    };

    for (cases) |entry| {
        const case, const code = entry;
        const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
        defer allocator.free(zfi_bytes);
        const line_count = std.mem.readInt(u64, zfi_bytes[table_offset..][0..8], .little);
        try std.testing.expect(line_count >= 2);
        const first_line = table_offset + @sizeOf(u64);
        const second_line = first_line + @sizeOf(main.index_format.SideTableLine);
        const last_line = first_line + (@as(usize, @intCast(line_count)) - 1) * @sizeOf(main.index_format.SideTableLine);

        switch (case) {
            .broken_base_prefix => std.mem.writeInt(u64, zfi_bytes[first_line..][0..8], 1, .little),
            .repeated_byte_offset => {
                const first_offset = std.mem.readInt(u64, zfi_bytes[first_line + 8 ..][0..8], .little);
                std.mem.writeInt(u64, zfi_bytes[second_line + 8 ..][0..8], first_offset, .little);
            },
            .overlapping_line_span => {
                const first_offset = std.mem.readInt(u64, zfi_bytes[first_line + 8 ..][0..8], .little);
                const first_bytes = std.mem.readInt(u32, zfi_bytes[first_line + 16 ..][0..4], .little);
                std.mem.writeInt(u64, zfi_bytes[second_line + 8 ..][0..8], first_offset + first_bytes - 1, .little);
            },
            .line_past_fasta => std.mem.writeInt(u32, zfi_bytes[first_line + 16 ..][0..4], fasta_data.len + 1, .little),
            .wrong_base_total => {
                const line_bases = std.mem.readInt(u32, zfi_bytes[last_line + 20 ..][0..4], .little);
                std.mem.writeInt(u32, zfi_bytes[last_line + 20 ..][0..4], line_bases + 1, .little);
            },
        }

        const paths = try writeFastaAndRawZfi(allocator, @tagName(case), fasta_data, zfi_bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

        try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, code);
    }
}

test "[failure] - [zfi loader]: rejects overlapping side tables" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">a\nAAA\nCCCC\nGG\n>b\nTT\nAAAA\nC\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
    try std.testing.expect(!index.records.items[1].isUniformWidth());

    index.records.items[1].seq_offset = index.records.items[0].seq_offset;
    index.records.items[1].seq_len = index.records.items[0].seq_len;
    try index.records.items[1].markNonUniform(index.records.items[0].sideTableOffset());

    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    defer allocator.free(zfi_bytes);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-side-overlap", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, .@"index.side_table_order_invalid");
    for ([_]main.index_format.LoadMode{ .catalog, .stats_scan }) |mode| {
        var failure: diagnostics.Failure = undefined;
        try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.side_table_order_invalid", &failure);
        try std.testing.expectEqual(@as(u64, 2), failure.primary().?.at.record.?);
    }
}

test "[integration] - [zfi loader]: accepts side tables and an embedded name blob" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    defer allocator.free(zfi_bytes);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-production-ok", fasta_data, zfi_bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var idx = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .on_access);
    defer idx.deinit();
    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    try std.testing.expect(idx.name_blob != null);
    try std.testing.expectEqual(@as(usize, 1), idx.records.len());
    try std.testing.expectEqualStrings("seq", idx.recordName(0).?);
}

test "[failure] - [index loader]: every allocation failure releases a partial load" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        loadAndDeinitForAllocationCheck,
        .{ "tests/data/simple.fasta", main.index_format.LoadMode.lookup_full_map },
    );

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var fasta = std.Io.Writer.Allocating.init(allocator);
    defer fasta.deinit();
    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    var last_offset: usize = 0;
    for (0..64) |i| {
        try fasta.writer.print(">r{d}\n", .{i});
        last_offset = fasta.written().len;
        try fasta.writer.writeAll("ACGT\n");
        try fai.writer.print("r{d}\t4\t{d}\t4\t5\n", .{ i, last_offset });
    }
    const paths = try writeFastaAndFai(allocator, "fai-allocation-failures", fasta.written(), fai.written());
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        loadAndDeinitForAllocationCheck,
        .{ paths.fasta_path, main.index_format.LoadMode.stats_scan },
    );
    var idx = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .stats_scan);
    defer idx.deinit();
    try std.testing.expectEqual(@as(usize, 64), idx.records.len());
    try std.testing.expectEqual(@as(u64, last_offset), (try recordAt(&idx, 63)).seq_offset);
    try std.testing.expectEqual(@as(u64, 4), (try recordAt(&idx, 63)).seq_len);
    try std.testing.expectEqual(@as(usize, 2), (try recordAt(&idx, 0)).name_len);
}

test "[property] - [zfi loader]: modes select the first duplicate" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">dup\nAAAA\n>dup\nCCCCCC\n";
    const fasta_path = try writeFastaArtifact(allocator, "duplicate-zfi", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    var full_map_idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer full_map_idx.deinit();
    const requested = [_][]const u8{"dup"};
    var positional_idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .{ .requested = &requested });
    defer positional_idx.deinit();

    try std.testing.expectEqual(full_map_idx.lookupName("dup"), positional_idx.lookupName("dup"));
    try std.testing.expectEqual(@as(?usize, 0), positional_idx.lookupName("dup"));
}

test "[integration] - [index loader]: gives a present ZFI authority over FAI" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGTACGT\n>seq2\nGGGG\n";
    const fasta_path = try writeFastaArtifact(allocator, "stale-zfi", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
    defer fai_file.close(io);
    var fai_buf: [65536]u8 = undefined;
    var fai_fw = fai_file.writer(io, &fai_buf);
    _ = try scanFaiData(allocator, fasta_data, &fai_fw.interface);
    try fai_fw.flush();

    {
        var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
        defer idx.deinit();
        try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    }

    const zfi_file = try std.Io.Dir.cwd().openFile(io, zfi_path, .{});
    defer zfi_file.close(io);
    try markFileStaleOneHourAgo(zfi_file);

    try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"index.older_than_source");
    try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta_path, "seq1" }, .@"index.older_than_source");
}

test "[failure] - [index loader]: rejects an invalid present ZFI despite valid FAI" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "invalid-zfi-authority", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, "seq\t4\t5\t4\t5\n");
    }

    const invalid_indexes = [_]struct { []const u8, diagnostics.Code }{
        .{ "", .@"index.truncated" },
        .{ "ZFI", .@"index.truncated" },
        .{ "ZFI\x03 incompatible", .@"index.unsupported_version" },
        .{ "not a zfi index", .@"index.truncated" },
    };
    for (invalid_indexes) |entry| {
        const bytes, const code = entry;
        {
            const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
            defer zfi_file.close(io);
            try std.Io.File.writeStreamingAll(zfi_file, io, bytes);
        }
        try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, code);
    }
}

test "[failure] - [index loader]: reports an inaccessible selected sidecar" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "inaccessible-sidecar", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{});
        defer zfi_file.close(io);
        try std.Io.File.writeStreamingAll(zfi_file, io, "not a zfi index");
    }
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, "seq\t4\t5\t4\t5\n");
    }

    const no_permissions: std.Io.File.Permissions = @enumFromInt(0);
    {
        const zfi_file = try std.Io.Dir.cwd().openFile(io, zfi_path, .{});
        defer zfi_file.close(io);
        const permissions = (try zfi_file.stat(io)).permissions;
        try zfi_file.setPermissions(io, no_permissions);
        defer zfi_file.setPermissions(io, permissions) catch {};
        try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"io.access_denied");
        try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", fasta_path, "seq" }, .@"io.access_denied");
    }

    try std.Io.Dir.cwd().deleteFile(io, zfi_path);
    {
        const fai_file = try std.Io.Dir.cwd().openFile(io, fai_path, .{});
        defer fai_file.close(io);
        const permissions = (try fai_file.stat(io)).permissions;
        try fai_file.setPermissions(io, no_permissions);
        defer fai_file.setPermissions(io, permissions) catch {};
        try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"io.access_denied");
    }
}

test "[property] - [fai loader]: modes select the first duplicate" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">dup\nAAAA\n>dup\nCCCCCC\n";
    const paths = try writeFastaAndFai(allocator, "duplicate-fai", fasta_data, "dup\t4\t5\t4\t5\ndup\t6\t15\t6\t7\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    var full_map_idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer full_map_idx.deinit();
    const requested = [_][]const u8{"dup"};
    var positional_idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .{ .requested = &requested });
    defer positional_idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, positional_idx.source);
    try std.testing.expectEqual(@as(usize, 1), positional_idx.records.len());
    try std.testing.expectEqualStrings("dup", positional_idx.recordName(0).?);
    try std.testing.expectEqual(full_map_idx.lookupName("dup"), positional_idx.lookupName("dup"));
    try std.testing.expectEqual(@as(?usize, 0), positional_idx.lookupName("dup"));
}

test "[integration] - [fai loader]: requested names retain first matches and check only the lines they keep" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(
        allocator,
        "positional-fai",
        ">wanted\nAAAA\n>wanted\nCCCCCC\n>colon:name\nGGGG\n>\nTT\n>other\nAC\n",
    );
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(
            fai_file,
            io,
            "wanted\t4\t8\t4\t5\nwanted\t6\t21\t6\t7\ncolon:name\t4\t40\t4\t5\n\t2\t47\t2\t3\nother\t2\t57\t2\t3\n",
        );
    }

    const requested = [_][]const u8{ "wanted", "missing", "colon:name", "" };
    {
        var idx = try support.loadIndex(
            std.testing.allocator,
            io,
            fasta_path,
            .{ .requested = requested[0..] },
        );
        defer idx.deinit();

        try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, idx.source);
        try std.testing.expect(idx.name_map != null);
        try std.testing.expectEqual(@as(usize, 3), idx.records.len());
        try std.testing.expectEqual(@as(usize, 3), idx.name_slices.len);
        try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("wanted"));
        try std.testing.expectEqual(@as(?usize, 1), idx.lookupName("colon:name"));
        try std.testing.expectEqual(@as(?usize, 2), idx.lookupName(""));
        try std.testing.expectEqual(@as(?usize, null), idx.lookupName("missing"));
        try std.testing.expectEqual(@as(u64, 4), (try recordAt(&idx, 0)).seq_len);
        try std.testing.expectEqualStrings("wanted", idx.recordName(0).?);
    }

    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(
            fai_file,
            io,
            "wanted\t4\t8\t4\t5\nwanted\t6\t21\t6\t7\ncolon:name\t4\t40\t4\t5\n\t2\t47\t2\t3\nother\t2\t57\t2\t1\n",
        );
    }
    {
        var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .{ .requested = requested[0..] });
        defer idx.deinit();
        try std.testing.expectEqual(@as(usize, 3), idx.records.len());
    }

    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(
            fai_file,
            io,
            "wanted\t4\t8\t4\t5\nwanted\t6\t21\t6\t7\ncolon:name\t4\t40\t4\t1\n\t2\t47\t2\t3\nother\t2\t57\t2\t3\n",
        );
    }
    try expectLoadFailure(std.testing.allocator, fasta_path, .{ .requested = requested[0..] }, .@"index.invalid_line_geometry");
}

test "[failure] - [fai loader]: any-name lookups check only the lines they use" {
    const allocator = std.testing.allocator;
    const paths = try writeFastaAndFai(allocator, "fai-lazy", ">a\nACGT\n>b\nGG\n>a\nT\n", "a\t4\t3\t4\t5\n\nb\t2\t11\t2\t1\na\t1\t16\t1\t2");
    defer allocator.free(paths.fasta_path);
    defer allocator.free(paths.fai_path);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
    var idx = try support.loadIndex(allocator, io, paths.fasta_path, .lookup_full_map);
    defer idx.deinit();
    try std.testing.expectEqual(@as(usize, 3), idx.records.len());
    try std.testing.expectEqual(@as(?usize, 0), idx.lookupName("a"));
    try std.testing.expectEqual(@as(?usize, 1), idx.lookupName("b"));
    try std.testing.expectEqualStrings("a", idx.recordName(2).?);
    try std.testing.expectEqual(@as(u64, 4), (try recordAt(&idx, 0)).seq_len);
    try std.testing.expectEqual(@as(u64, 1), (try recordAt(&idx, 2)).seq_len);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectError(error.InvalidIndex, idx.checkedRecord(1, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(u64, 3), primary.at.location.line.?);
    try std.testing.expectEqual(@as(u64, 2), primary.at.record.?);
}

test "[unit] - [fai loader]: full lookup borrows mapped names" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "borrowed-fai-name", ">seq\nACGT\n", "seq\t4\t5\t4\t5\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer idx.deinit();
    const mapped_name = (try recordAt(&idx, 0)).getName(idx.fai_data.?);
    var names = idx.name_map.?.iterator();
    const map_entry = names.next().?;

    try std.testing.expectEqual(@as(usize, 0), idx.name_slices.len);
    try std.testing.expectEqual(mapped_name.ptr, idx.recordName(0).?.ptr);
    try std.testing.expectEqual(mapped_name.ptr, map_entry.key_ptr.*.ptr);
}

test "[integration] - [fai loader]: a requested-names load keeps the FASTA descriptor" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "fai-source-descriptor", ">seq\nACGT\n", "seq\t4\t5\t4\t5\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    const requested = [_][]const u8{"seq"};
    var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .{ .requested = &requested });
    defer idx.deinit();

    var first_byte: [1]u8 = undefined;
    try std.testing.expectEqual(
        @as(usize, 1),
        try std.Io.File.readPositionalAll(idx.fasta_file, io, &first_byte, 0),
    );
    try std.testing.expectEqual(@as(u8, '>'), first_byte[0]);
}

test "[integration] - [index loader]: keeps FASTA readable without mapping it" {
    var zfi_idx = try support.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta", .on_access);
    defer zfi_idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, zfi_idx.source);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "fai-get-source", ">seq\nACGT\n", "seq\t4\t5\t4\t5\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    var fai_idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer fai_idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.fai, fai_idx.source);

    var first_byte: [1]u8 = undefined;
    try std.testing.expectEqual(
        @as(usize, 1),
        try std.Io.File.readPositionalAll(fai_idx.fasta_file, io, &first_byte, 0),
    );
    try std.testing.expectEqual(@as(u8, '>'), first_byte[0]);
}

test "[integration] - [index loader]: keeps ZFI FASTA readable without mapping it" {
    var idx = try support.loadIndex(std.testing.allocator, io, "tests/data/simple.fasta", .stats_scan);
    defer idx.deinit();

    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    try std.testing.expect(idx.zfi_map != null);
    try std.testing.expectEqualStrings("seq1", idx.recordName(0).?);

    var first_byte: [1]u8 = undefined;
    try std.testing.expectEqual(
        @as(usize, 1),
        try std.Io.File.readPositionalAll(idx.fasta_file, io, &first_byte, 0),
    );
    try std.testing.expectEqual(@as(u8, '>'), first_byte[0]);
}

test "[failure] - [fai loader]: rejects malformed numeric fields and impossible geometry" {
    const Case = struct {
        stem: []const u8,
        fai: []const u8,
        code: diagnostics.Code,
    };
    const cases = [_]Case{
        .{ .stem = "bad-fai-zero-length", .fai = "seq1\t0\t6\t4\t5\n", .code = .@"index.zero_sequence_length" },
        .{ .stem = "bad-fai-zero-bases", .fai = "seq1\t4\t6\t0\t5\n", .code = .@"index.invalid_line_geometry" },
        .{ .stem = "bad-fai-zero-bytes", .fai = "seq1\t4\t6\t4\t0\n", .code = .@"index.invalid_line_geometry" },
        .{ .stem = "bad-fai-line-bytes", .fai = "seq1\t4\t6\t5\t4\n", .code = .@"index.invalid_line_geometry" },
        .{ .stem = "bad-fai-offset", .fai = "seq1\t4\t11\t4\t5\n", .code = .@"index.sequence_offset_out_of_bounds" },
        .{ .stem = "bad-fai-span", .fai = "seq1\t100\t6\t4\t5\n", .code = .@"index.sequence_span_out_of_bounds" },
        .{ .stem = "bad-fai-max-length", .fai = "seq1\t18446744073709551615\t6\t4\t5\n", .code = .@"index.sequence_span_out_of_bounds" },
        .{ .stem = "bad-fai-u64-overflow", .fai = "seq1\t99999999999999999999\t6\t4\t5\n", .code = .@"index.fai.integer_out_of_range" },
        .{ .stem = "bad-fai-u32-overflow", .fai = "seq1\t4\t6\t4294967296\t5\n", .code = .@"index.fai.integer_out_of_range" },
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (cases) |case| {
        const paths = try writeFastaAndFai(allocator, case.stem, ">seq1\nACGT\n", case.fai);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        try expectLoadFailure(std.testing.allocator, paths.fasta_path, .catalog, case.code);
    }
}

test "[property] - [fai loader]: rejects malformed fields in every mode" {
    const Case = struct {
        stem: []const u8,
        fai: []const u8,
        code: diagnostics.Code,
    };
    const cases = [_]Case{
        .{ .stem = "empty-fai", .fai = "", .code = .@"index.no_records" },
        .{ .stem = "blank-fai", .fai = "\n\n", .code = .@"index.no_records" },
        .{ .stem = "missing-fields", .fai = "seq1\t4\t6\t4\n", .code = .@"index.fai.field_missing" },
        .{ .stem = "empty-field", .fai = "seq1\t\t6\t4\t5\n", .code = .@"index.fai.field_missing" },
        .{ .stem = "non-decimal-field", .fai = "seq1\t4x\t6\t4\t5\n", .code = .@"index.fai.invalid_integer" },
        .{ .stem = "signed-field", .fai = "seq1\t-4\t6\t4\t5\n", .code = .@"index.fai.invalid_integer" },
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (cases) |case| {
        const paths = try writeFastaAndFai(allocator, case.stem, ">seq1\nACGT\n", case.fai);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        for ([_]main.index_format.LoadMode{ .catalog, .{ .requested = &.{"seq1"} }, .stats_scan }) |mode| {
            try expectLoadFailure(std.testing.allocator, paths.fasta_path, mode, case.code);
        }
    }
}

test "[failure] - [fai loader]: a requested-names load rejects impossible geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "bad-fai-positional", ">seq1\nACGT\n", "seq1\t4\t6\t0\t5\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    try expectLoadFailure(std.testing.allocator, fasta_path, .{ .requested = &.{"seq1"} }, .@"index.invalid_line_geometry");
}

test "[failure] - [fai loader]: stats mode rejects impossible geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const paths = try writeFastaAndFai(allocator, "bad-fai-stats-scan", ">seq1\nACGT\n", "seq1\t4\t6\t4\t3\n");
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    try expectLoadFailure(std.testing.allocator, fasta_path, .stats_scan, .@"index.invalid_line_geometry");
}

test "[failure] - [fai loader]: rejects names above the u16 limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fai_buf: std.ArrayList(u8) = .empty;
    defer fai_buf.deinit(allocator);
    try fai_buf.appendNTimes(allocator, 'N', 65_600);
    try fai_buf.appendSlice(allocator, "\t4\t6\t4\t5\n");
    const paths = try writeFastaAndFai(allocator, "bad-fai-long-name", ">seq1\nACGT\n", fai_buf.items);
    const fasta_path = paths.fasta_path;
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    for ([_]main.index_format.LoadMode{ .catalog, .{ .requested = &.{"seq1"} }, .stats_scan }) |mode| {
        try expectLoadFailure(std.testing.allocator, fasta_path, mode, .@"index.identifier_too_long");
    }
}

test "[edge] - [fai loader]: streamed modes accept the maximum sequence name" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try uniqueArtifactPath(allocator, "fai-max-name", "fa");
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};

    const name = try allocator.alloc(u8, std.math.maxInt(u16));
    @memset(name, 'n');

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, ">");
        try std.Io.File.writeStreamingAll(fasta_file, io, name);
        try std.Io.File.writeStreamingAll(fasta_file, io, "\nA\n");
    }
    {
        const fai_file = try std.Io.Dir.cwd().createFile(io, fai_path, .{ .truncate = true });
        defer fai_file.close(io);
        try std.Io.File.writeStreamingAll(fai_file, io, name);
        try std.Io.File.writeStreamingAll(fai_file, io, "\t1\t65537\t1\t2\n");
    }

    {
        var idx = try support.loadIndex(
            std.testing.allocator,
            io,
            fasta_path,
            .{ .requested = &.{name} },
        );
        defer idx.deinit();

        try std.testing.expectEqual(@as(usize, 1), idx.records.len());
        try std.testing.expectEqualStrings(name, idx.recordName(0).?);
    }
    {
        var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .stats_scan);
        defer idx.deinit();

        try std.testing.expectEqual(@as(usize, 1), idx.records.len());
        try std.testing.expectEqual(name.len, (try recordAt(&idx, 0)).name_len);
    }
}

test "[property] - [fai loader]: modes agree when FASTA or index omits the terminal newline" {
    const Case = struct {
        stem: []const u8,
        fasta: []const u8,
        fai: []const u8,
    };
    const cases = [_]Case{
        .{
            .stem = "fai-no-term-nl",
            .fasta = ">s1\nACGT\n>s2\nGGGGGG\n>s3\nCCCCC\n",
            .fai = "s1\t4\t4\t4\t5\ns2\t6\t13\t6\t7\ns3\t5\t24\t5\t6",
        },
        .{
            .stem = "fasta-no-term-nl",
            .fasta = ">seq1\nACGT",
            .fai = "seq1\t4\t6\t4\t5\n",
        },
        .{
            .stem = "both-no-term-nl",
            .fasta = ">seq1\nAAAA\nBBBB",
            .fai = "seq1\t8\t6\t4\t5",
        },
    };
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (cases) |case| {
        const paths = try writeFastaAndFai(allocator, case.stem, case.fasta, case.fai);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        try expectFaiLoaderModesAgree(allocator, paths.fasta_path);
    }
}

test "[edge] - [zfi loader]: accepts FASTA without a terminal newline" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT";
    const fasta_path = try writeFastaArtifact(allocator, "zfi-fasta-no-term-nl", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    var idx = try support.loadIndex(std.testing.allocator, io, fasta_path, .on_access);
    defer idx.deinit();
    try std.testing.expectEqual(main.index_format.LoadedIndex.IndexSource.zfi, idx.source);
    try std.testing.expectEqual(@as(u64, 4), (try recordAt(&idx, 0)).seq_len);
}

test "[failure] - [fai scanner]: rejects non-uniform sequence layout without output" {
    const data =
        \\>mixed
        \\AAAA
        \\BBBBBB
        \\CC
        \\
    ;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var out = std.Io.Writer.Allocating.init(allocator);
    defer out.deinit();
    try std.testing.expectError(
        error.InvalidIndex,
        scanFaiData(allocator, data, &out.writer),
    );
    try std.testing.expectEqual(@as(usize, 0), out.written().len);
}

test "[integration] - [fai output]: uses a separate spool and cleans it after success" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const source_dir_path = try uniqueArtifactDirPath(allocator, "fai-read-only-source");
    defer std.Io.Dir.cwd().deleteTree(io, source_dir_path) catch {};
    const spool_dir_path = try uniqueArtifactDirPath(allocator, "fai-spool-success");
    defer std.Io.Dir.cwd().deleteTree(io, spool_dir_path) catch {};
    const fasta_path = try std.fmt.allocPrint(allocator, "{s}/records.fa", .{source_dir_path});

    var expected = std.Io.Writer.Allocating.init(allocator);
    defer expected.deinit();
    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{});
        defer fasta_file.close(io);
        var file_buf: [4096]u8 = undefined;
        var fasta_writer = fasta_file.writer(io, &file_buf);
        for (0..2048) |record_index| {
            const offset = fasta_writer.logicalPos();
            try fasta_writer.interface.print(">r{d}\nA\n", .{record_index});
            try expected.writer.print("r{d}\t1\t{d}\t1\t2\n", .{
                record_index,
                offset + 3 + std.fmt.count("{d}", .{record_index}),
            });
        }
        try fasta_writer.interface.flush();
    }

    var source_dir = try std.Io.Dir.cwd().openDir(io, source_dir_path, .{ .iterate = true });
    defer source_dir.close(io);
    const source_permissions = (try source_dir.stat(io)).permissions;
    try source_dir.setPermissions(io, source_permissions.setReadOnly(true));
    defer source_dir.setPermissions(io, source_permissions) catch {};

    var env = try std.testing.environ.createMap(allocator);
    defer env.deinit();
    try env.put("TMPDIR", spool_dir_path);
    const missing_fallback = try std.fmt.allocPrint(allocator, "{s}/missing", .{spool_dir_path});
    try env.put("TEMP", missing_fallback);
    try env.put("TMP", missing_fallback);

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "index", "--emit-fai", fasta_path },
        .environ_map = &env,
        .stdout_limit = .limited(1024 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| if (code != 0) {
            std.debug.print("emit-fai subprocess exited {d}: {s}\n", .{ code, result.stderr });
            return error.ChildProcessFailed;
        },
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqualSlices(u8, expected.written(), result.stdout);
    try std.testing.expect(result.stdout.len > (try std.Io.Dir.cwd().statFile(io, fasta_path, .{})).size);
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    try expectDirectoryEmpty(spool_dir_path);

    var sibling_buf: [4096]u8 = undefined;
    const sibling = try std.fmt.bufPrint(&sibling_buf, "{s}.fai.stdout.tmp", .{fasta_path});
    try std.testing.expectError(error.FileNotFound, std.Io.Dir.cwd().access(io, sibling, .{}));
}

test "[failure] - [fai output]: leaves output and spool empty" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(
        allocator,
        "fai-nonuniform-spool",
        ">mixed\nAAAA\nBBBBBB\nCC\n",
    );
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const spool_dir_path = try uniqueArtifactDirPath(allocator, "fai-spool-failure");
    defer std.Io.Dir.cwd().deleteTree(io, spool_dir_path) catch {};
    var env = try std.testing.environ.createMap(allocator);
    defer env.deinit();
    try env.put("TMPDIR", "zig-cache/test-artifacts/definitely-missing-spool");
    try env.put("TEMP", spool_dir_path);
    try env.put("TMP", "");

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{
        .argv = &.{ ZFASTA_BIN, "index", "--error-format=json", "--emit-fai", fasta_path },
        .environ_map = &env,
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);

    switch (result.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try support.expectDiagnostic(allocator, result.stderr, .@"index.nonuniform_fai");
    try expectDirectoryEmpty(spool_dir_path);
}

test "[failure] - [zfi loader]: rejects a source size change" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try writeFastaArtifact(allocator, "stale-size", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{ .truncate = true });
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, ">seq1\nACGTT\n");
    }
    try keepIndexNewer(zfi_path, fasta_path);

    try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"index.source_size_mismatch");
}

test "[failure] - [zfi loader]: rejects a same-size source replacement" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_a = ">seq1\nACGT\n";
    const fasta_b = ">seq1\nTTTT\n";
    const fasta_path = try writeFastaArtifact(allocator, "stale-same-size", fasta_a);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    try std.testing.expectEqual(fasta_a.len, fasta_b.len);
    const mtime_a = try statMtimeNs(fasta_path);
    var index = try main.indexer.scanZfiData(allocator, fasta_a);
    defer index.deinit(allocator);
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_a.len, mtime_a);

    {
        const fasta_file = try std.Io.Dir.cwd().createFile(io, fasta_path, .{ .truncate = true });
        defer fasta_file.close(io);
        try std.Io.File.writeStreamingAll(fasta_file, io, fasta_b);
        const bumped = std.Io.Timestamp.fromNanoseconds(@as(i96, @intCast(mtime_a + std.time.ns_per_s)));
        try fasta_file.setTimestamps(io, .{
            .modify_timestamp = .{ .new = bumped },
        });
    }
    try keepIndexNewer(zfi_path, fasta_path);

    try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"index.source_mtime_mismatch");
}

test "[integration] - [zfi loader]: names that end like the source identity load exactly" {
    const allocator = std.testing.allocator;
    const name = "nZFIDabcdefgh";
    for ([_]usize{ 0, 1800 }) |prefix_count| {
        for ([_][]const u8{ ">nZFIDabcdefgh\nACGT\n", ">nZFIDabcdefgh\nACG\nTACGT\nAC\n" }) |sample| {
            var input: std.ArrayList(u8) = .empty;
            defer input.deinit(allocator);
            for (0..prefix_count) |i| {
                var header: [64]u8 = undefined;
                try input.appendSlice(allocator, try std.fmt.bufPrint(&header, ">prefix{d}\nACGT\n", .{i}));
            }
            try input.appendSlice(allocator, sample);
            const fasta = input.items;
            var catalog = try main.indexer.scanZfiData(allocator, fasta);
            defer catalog.deinit(allocator);
            {
                const path = try writeFastaArtifact(allocator, "zfi-name-identity-collision", fasta);
                defer allocator.free(path);
                defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
                const bytes = try zfiBytes(allocator, &catalog, fasta.len, try statMtimeNs(path));
                defer allocator.free(bytes);
                const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
                defer allocator.free(zfi_path);
                defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
                {
                    const file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{});
                    defer file.close(io);
                    try std.Io.File.writeStreamingAll(file, io, bytes);
                }
                const modes = [_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{name} }, .stats_scan };
                for (modes) |mode| {
                    var failure: diagnostics.Failure = undefined;
                    failure.reset(.get);
                    var loaded = try main.index_format.loadIndex(allocator, io, path, .automatic, mode, &failure);
                    defer loaded.deinit();
                    const last = loaded.records.len() - 1;
                    try std.testing.expectEqual(@as(?usize, last), try loaded.findName(name, &failure));
                    const rec = try loaded.checkedRecord(last, &failure);
                    try std.testing.expectEqual(catalog.records.items[last].seq_len, rec.seq_len);
                    var output = std.Io.Writer.Allocating.init(allocator);
                    defer output.deinit();
                    try main.getter.extractRegion(&loaded, name, &output.writer, &failure);
                    const expected = if (rec.seq_len == 4) ">nZFIDabcdefgh\nACGT\n" else ">nZFIDabcdefgh\nACGTACGTAC\n";
                    try std.testing.expectEqualStrings(expected, output.written());
                }
            }
        }
    }
}

test "[failure] - [zfi loader]: rejects an index without source identity" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta_data = ">seq1\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, 0);
    const footer_len = main.index_format.ZFI_NAME_FOOTER_BYTES;
    const id_off = zfi_bytes.len - footer_len - main.index_format.ZFI_SOURCE_ID_BYTES;
    var without_identity: std.ArrayList(u8) = .empty;
    try without_identity.appendSlice(allocator, zfi_bytes[0..id_off]);
    try without_identity.appendSlice(allocator, zfi_bytes[zfi_bytes.len - footer_len ..]);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-no-identity", fasta_data, without_identity.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"seq1"} } }) |mode| {
        var failure: diagnostics.Failure = undefined;
        try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.name_footer_missing", &failure);
    }
}

test "[failure] - [zfi loader]: inflated blob length cannot hide source identity" {
    const allocator = std.testing.allocator;
    const fasta = ">nZFIDabcdefgh\nACGT\n";
    var catalog = try main.indexer.scanZfiData(allocator, fasta);
    defer catalog.deinit(allocator);
    const bytes = try zfiBytes(allocator, &catalog, fasta.len, 0);
    defer allocator.free(bytes);
    std.mem.writeInt(u64, bytes[bytes.len - 8 ..][0..8], catalog.name_blob.items.len + main.index_format.ZFI_SOURCE_ID_BYTES, .little);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-inflated-identity-blob", fasta, bytes);
    defer allocator.free(paths.fasta_path);
    defer allocator.free(paths.zfi_path);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"nZFIDabcdefgh"} } }) |mode| {
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        try std.testing.expectError(error.InvalidIndex, main.index_format.loadIndex(allocator, io, paths.fasta_path, .automatic, mode, &failure));
        try std.testing.expectEqual(diagnostics.Code.@"index.name_blob_out_of_bounds", failure.primary().?.code());
    }
}

test "[failure] - [zfi loader]: rejects mismatched embedded source time" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try writeFastaArtifact(allocator, "zfi-bad-mtime", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    const real_mtime = try statMtimeNs(fasta_path);
    const wrong_mtime = real_mtime +% 1;

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    const zfi_bytes = try zfiBytes(allocator, &index, fasta_data.len, wrong_mtime);
    defer allocator.free(zfi_bytes);

    {
        const zfi_file = try std.Io.Dir.cwd().createFile(io, zfi_path, .{ .truncate = true });
        defer zfi_file.close(io);
        try std.Io.File.writeStreamingAll(zfi_file, io, zfi_bytes);
    }

    try expectLoadFailure(std.testing.allocator, fasta_path, .on_access, .@"index.source_mtime_mismatch");
}

test "[failure] - [zfi loader]: rejects truncation inside a record" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var zfi_buf: std.ArrayList(u8) = .empty;
    defer zfi_buf.deinit(allocator);
    try appendZfiHeader(&zfi_buf, allocator, 1, 11);
    try zfi_buf.appendSlice(allocator, &[_]u8{0} ** 20);

    const paths = try writeFastaAndRawZfi(allocator, "zfi-trunc-record", ">seq1\nACGT\n", zfi_buf.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    try expectLoadFailure(std.testing.allocator, paths.fasta_path, .on_access, .@"index.truncated");
}

test "[failure] - [zfi loader]: rejects zero line geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_data = ">seq1\nACGT\n";
    const fasta_path = try writeFastaArtifact(allocator, "zfi-zero-geom", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    index.records.items[0].line_bases = 0;
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));
    try expectLoadFailure(std.testing.allocator, fasta_path, .catalog, .@"index.invalid_line_geometry");
}

test "[edge] - [zfi scanner]: keeps a trailing-blank record uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, ">seq\nACGT\n\n");
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
}

test "[edge] - [fasta scanner]: ignores a trailing blank before the next header" {
    const data = ">a\nAAAA\n\n>b\nCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expect(index.records.items[1].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
    const second_name = index.name_blob.items[index.records.items[1].name_offset..][0..index.records.items[1].name_len];
    try std.testing.expectEqualStrings("b", second_name);

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    const count = try scanFaiData(allocator, data, &fai.writer);
    try std.testing.expectEqual(@as(u32, 2), count);
    try std.testing.expectEqualStrings("a\t4\t3\t4\t5\nb\t4\t12\t4\t5\n", fai.written());
}

test "[edge] - [zfi scanner]: marks an interior blank non-uniform" {
    const data = ">seq\nAAAA\n\nCCCC\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
}

test "[edge] - [fasta scanner]: rejects an interior blank after a uniform stride prefix" {
    const data = ">seq\nAAAA\nAAAA\nAAAA\n\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(
        @as(u64, 4),
        std.mem.readInt(u64, index.side_tables.items[0..@sizeOf(u64)], .little),
    );

    var fai = std.Io.Writer.Allocating.init(allocator);
    defer fai.deinit();
    try std.testing.expectError(
        error.InvalidIndex,
        scanFaiData(allocator, data, &fai.writer),
    );
}

test "[edge] - [zfi scanner]: keeps many full lines and a short final line uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">chrom\n");
    var line: [61]u8 = undefined;
    @memset(line[0..60], 'A');
    line[60] = '\n';
    var i: usize = 0;
    while (i < 4000) : (i += 1) {
        try fasta.appendSlice(allocator, line[0..]);
    }
    try fasta.appendSlice(allocator, "ACGT\n");

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), index.records.items.len);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
    try std.testing.expectEqual(@as(u64, 4000 * 60 + 4), index.records.items[0].seq_len);
}

test "[edge] - [zfi scanner]: creates a side table for an interior short line" {
    const data = ">seq\nAAAA\nAAAA\nAA\nAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
    try std.testing.expect(index.side_tables.items.len > 0);
}

test "[edge] - [zfi scanner]: marks a long final line non-uniform" {
    const data = ">seq\nAAAA\nAAAA\nAAAAAA\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(!index.records.items[0].isUniformWidth());
    try std.testing.expect(index.side_tables.items.len > 0);
}

test "[edge] - [zfi scanner]: keeps a short final line without newline uniform" {
    const data = ">seq\nAAAA\nAAAA\nAC";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
}

test "[edge] - [zfi scanner]: keeps a short final line across the read boundary uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq\n");
    const full_line = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    while (fasta.items.len + full_line.len < INDEX_READ_BUFFER_BYTES + 200) {
        try fasta.appendSlice(allocator, full_line);
    }
    try fasta.appendSlice(allocator, "ACGT\n");

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
}

test "[edge] - [zfi scanner]: keeps a CRLF body with a short final line uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">chrom\r\n");
    var line: [62]u8 = undefined;
    @memset(line[0..60], 'A');
    line[60] = '\r';
    line[61] = '\n';
    var i: usize = 0;
    while (i < 200) : (i += 1) {
        try fasta.appendSlice(allocator, line[0..]);
    }
    try fasta.appendSlice(allocator, "ACGT\r\n");

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expect(index.records.items[0].isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), index.side_tables.items.len);
    try std.testing.expectEqual(@as(u32, 62), index.records.items[0].line_bytes);
}

test "[property] - [fasta scanner]: stride blocks fall back at validation boundaries" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const full_line = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    var read_storage: [64 * 1024]u8 = undefined;
    for ([_]usize{ 255, 256, 257 }) |break_line| {
        var fasta: std.ArrayList(u8) = .empty;
        defer fasta.deinit(allocator);
        try fasta.appendSlice(allocator, ">seq\n");
        for (0..520) |line_index| {
            try fasta.appendSlice(
                allocator,
                if (line_index == break_line) "AAAA\n" else full_line,
            );
        }

        {
            try expectReaderMatchesProduction(
                allocator,
                fasta.items,
                &read_storage,
                "stride-block-boundary",
            );
        }
    }
}

test "[property] - [fasta scanner]: CRLF stride blocks preserve both formats" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">seq\r\n");
    const full_line = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\r\n";
    for (0..600) |_| {
        try fasta.appendSlice(allocator, full_line);
    }
    try fasta.appendSlice(allocator, "ACGT\r\n");

    var read_storage: [64 * 1024]u8 = undefined;
    {
        try expectReaderMatchesProduction(
            allocator,
            fasta.items,
            &read_storage,
            "crlf-stride-blocks",
        );
    }
}

test "[property] - [zfi scanner]: keeps multiple uniform bodies uniform" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    const full = "AAAA\n";
    var r: usize = 0;
    while (r < 3) : (r += 1) {
        var hdr_buf: [16]u8 = undefined;
        const hdr = try std.fmt.bufPrint(&hdr_buf, ">r{d}\n", .{r});
        try fasta.appendSlice(allocator, hdr);
        var i: usize = 0;
        while (i < 50) : (i += 1) {
            try fasta.appendSlice(allocator, full);
        }
        try fasta.appendSlice(allocator, "AC\n");
    }

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 3), index.records.items.len);
    for (index.records.items) |rec| {
        try std.testing.expect(rec.isUniformWidth());
    }
}

test "[regression] - [zfi scanner]: preserves a header spanning the prior stride width" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var fasta: std.ArrayList(u8) = .empty;
    defer fasta.deinit(allocator);
    try fasta.appendSlice(allocator, ">a\n");
    const full = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    var i: usize = 0;
    while (i < 5) : (i += 1) {
        try fasta.appendSlice(allocator, full);
    }
    try fasta.appendSlice(allocator, "NN\n");
    const hdr = ">b xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    try std.testing.expectEqual(@as(usize, 58), hdr.len);
    try fasta.appendSlice(allocator, hdr);
    try fasta.appendSlice(allocator, full);
    try fasta.appendSlice(allocator, "AC\n");

    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    const n0 = index.name_blob.items[index.records.items[0].name_offset..][0..index.records.items[0].name_len];
    const n1 = index.name_blob.items[index.records.items[1].name_offset..][0..index.records.items[1].name_len];
    try std.testing.expectEqualStrings("a", n0);
    try std.testing.expectEqualStrings("b", n1);
}

// --- Stable gates fixtures ---
test "[integration] - [fasta scanner]: duplicate names keep every occurrence" {
    try requireFixturePath("tests/data/gates/duplicates.fasta");
    const data = try readTestFile(std.testing.allocator, "tests/data/gates/duplicates.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 3), index.records.items.len);
    for (index.records.items, [_][]const u8{ "dup", "other", "dup" }) |record, name| {
        try std.testing.expectEqualStrings(name, index.name_blob.items[record.name_offset..][0..record.name_len]);
        try std.testing.expectEqual(@as(u64, 4), record.seq_len);
    }
}

test "[integration] - [fasta scanner]: enforces the header length boundary on the gate fixtures" {
    try requireFixturePath("tests/data/gates/long_header.fasta");
    try requireFixturePath("tests/data/gates/long_header_reject.fasta");

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const ok_data = try readTestFile(allocator, "tests/data/gates/long_header.fasta");
    defer allocator.free(ok_data);
    var ok_index = try main.indexer.scanZfiData(allocator, ok_data);
    defer ok_index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 1), ok_index.records.items.len);
    try std.testing.expectEqual(@as(u16, 1025), ok_index.records.items[0].name_len);

    const reject_data = try readTestFile(allocator, "tests/data/gates/long_header_reject.fasta");
    defer allocator.free(reject_data);
    try std.testing.expectEqual(@as(u8, '>'), reject_data[0]);
    try std.testing.expectEqual(null, std.mem.findScalar(u8, reject_data, '\r'));
    try std.testing.expectEqual(main.indexer.MAX_INDEX_NAME_LEN + 1, std.mem.findScalar(u8, reject_data, '\n').? - 1);
    try std.testing.expectError(error.InvalidIndex, main.indexer.scanZfiData(allocator, reject_data));
}

test "[integration] - [fasta scanner]: preserves invalid UTF-8 header bytes" {
    try requireFixturePath("tests/data/gates/invalid_utf8_header.fasta");
    const data = try readTestFile(std.testing.allocator, "tests/data/gates/invalid_utf8_header.fasta");
    defer std.testing.allocator.free(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var index = try main.indexer.scanZfiData(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), index.records.items.len);
    for (index.records.items) |record| {
        try std.testing.expectEqualSlices(u8, "bad\xffname", index.name_blob.items[record.name_offset..][0..record.name_len]);
    }
}

test "[integration] - [fai loader]: rejects the invalid-geometry gate fixtures" {
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta");
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta.fai");
    try requireFixturePath("tests/data/gates/seq1_large_offset.fasta");
    try requireFixturePath("tests/data/gates/seq1_large_offset.fasta.fai");

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectCorruptFaiFixture(allocator, "tests/data/gates/seq1_zero_geometry", .@"index.invalid_line_geometry");
    try expectCorruptFaiFixture(allocator, "tests/data/gates/seq1_large_offset", .@"index.sequence_offset_out_of_bounds");
}

test "[integration] - [zfi loader]: rejects zero geometry from the gate fixture" {
    try requireFixturePath("tests/data/gates/seq1_zero_geometry.fasta");
    const fasta_data = try readTestFile(std.testing.allocator, "tests/data/gates/seq1_zero_geometry.fasta");
    defer std.testing.allocator.free(fasta_data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const fasta_path = try writeFastaArtifact(allocator, "gates-zfi-zero", fasta_data);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};

    var index = try main.indexer.scanZfiData(allocator, fasta_data);
    defer index.deinit(allocator);
    index.records.items[0].line_bases = 0;
    try main.indexer.writeZfiIndexFile(io, zfi_path, &index, fasta_data.len, try statMtimeNs(fasta_path));
    try expectLoadFailure(std.testing.allocator, fasta_path, .catalog, .@"index.invalid_line_geometry");
}

test "[failure] - [index loader]: a directory sidecar fails as unmappable in every mode, empty or not" {
    const allocator = std.testing.allocator;
    const fasta_path = try writeFastaArtifact(allocator, "directory-sidecar", ">a\nACGT\n");
    defer allocator.free(fasta_path);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    for ([_][]const u8{ ".zfi", ".fai" }) |suffix| {
        const index_path = try std.fmt.allocPrint(allocator, "{s}{s}", .{ fasta_path, suffix });
        defer allocator.free(index_path);
        try std.Io.Dir.cwd().createDirPath(io, index_path);
        defer std.Io.Dir.cwd().deleteDir(io, index_path) catch {};
        const child_path = try std.fmt.allocPrint(allocator, "{s}/entry", .{index_path});
        defer allocator.free(child_path);
        for ([_]bool{ false, true }) |non_empty| {
            if (non_empty) {
                const child = try std.Io.Dir.cwd().createFile(io, child_path, .{});
                child.close(io);
            }
            defer if (non_empty) std.Io.Dir.cwd().deleteFile(io, child_path) catch {};
            for ([_]main.index_format.LoadMode{ .lookup_full_map, .catalog, .stats_scan, .{ .requested = &.{"a"} }, .{ .requested = &.{ "a", "b" } }, .on_access }) |mode| {
                var failure: diagnostics.Failure = undefined;
                failure.reset(.get);
                try std.testing.expectError(error.Io, main.index_format.loadIndex(allocator, io, fasta_path, .automatic, mode, &failure));
                try std.testing.expectEqual(diagnostics.Code.@"io.map_failed", failure.primary().?.code());
                try std.testing.expectEqual(null, failure.primary().?.cause);
            }
        }
    }
}

test "[integration] - [zfi loader]: lookups find first matches through the name table" {
    const allocator = std.testing.allocator;
    const fasta = ">other\nAC\nGTA\n>wanted\nAAA\nCCCC\nGG\n>wanted\nTTTT\n>uniform\nACTG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    defer allocator.free(bytes);
    const paths = try writeFastaAndRawZfi(allocator, "matched-zfi", fasta, bytes);
    defer allocator.free(paths.fasta_path);
    defer allocator.free(paths.zfi_path);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var full = try support.loadIndex(allocator, io, paths.fasta_path, .catalog);
    defer full.deinit();
    try std.testing.expect(full.on_demand == null);
    const modes = [_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"wanted"} }, .{ .requested = &.{} }, .{ .requested = &.{ "uniform", "wanted", "missing" } } };
    for (modes) |mode| {
        var loaded = try support.loadIndex(allocator, io, paths.fasta_path, mode);
        defer loaded.deinit();
        try std.testing.expect(loaded.on_demand != null and loaded.zfi_map == null);
        try std.testing.expectEqual(@as(usize, 4), loaded.records.len());
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        for ([_][]const u8{ "other", "wanted", "uniform", "missing", "" }, [_]?usize{ 0, 1, 3, null, null }) |name, expected| {
            try std.testing.expectEqual(expected, loaded.lookupName(name));
            try std.testing.expectEqual(expected, try loaded.findName(name, &failure));
            try std.testing.expectEqual(expected, full.lookupName(name));
        }
        for (0..4) |i| {
            const rec = try loaded.checkedRecord(i, &failure);
            try std.testing.expectEqual(try recordAt(&full, i), rec);
            try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(full.sideTableLines(rec)), std.mem.sliceAsBytes(loaded.sideTableLines(rec)));
        }
    }
    try std.testing.checkAllAllocationFailures(allocator, loadAndDeinitForAllocationCheck, .{ paths.fasta_path, main.index_format.LoadMode{ .requested = &.{ "wanted", "uniform" } } });
    try std.testing.checkAllAllocationFailures(allocator, loadAndDeinitForAllocationCheck, .{ paths.fasta_path, main.index_format.LoadMode.catalog });

    index.records.items[0].seq_len = 0;
    try main.indexer.writeZfiIndexFile(io, paths.zfi_path, &index, fasta.len, try statMtimeNs(paths.fasta_path));
    try expectLoadFailure(allocator, paths.fasta_path, .catalog, .@"index.zero_sequence_length");
    var lazy = try support.loadIndex(allocator, io, paths.fasta_path, .{ .requested = &.{"uniform"} });
    defer lazy.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectEqual(@as(?usize, 3), try lazy.findName("uniform", &failure));
    _ = try lazy.checkedRecord(3, &failure);
    try std.testing.expectError(error.InvalidIndex, lazy.checkedRecord(0, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"index.zero_sequence_length", failure.primary().?.code());
    try std.testing.expectEqual(@as(?u64, 1), failure.primary().?.at.record);
}

test "[edge] - [fai loader]: accepts long decimal fields and ignored extra columns" {
    const allocator = std.testing.allocator;
    var text = std.Io.Writer.Allocating.init(allocator);
    defer text.deinit();
    try text.writer.writeAll("a\t");
    try text.writer.splatByteAll('0', 70_000);
    try text.writer.writeAll("4\t3\t4\t5\t");
    try text.writer.splatByteAll('x', 70_000);
    try text.writer.writeAll("\nb\t2\t11\t2\t3\n");
    const paths = try writeFastaAndFai(allocator, "matched-long-fai", ">a\nACGT\n>b\nTT\n", text.written());
    defer allocator.free(paths.fasta_path);
    defer allocator.free(paths.fai_path);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};

    var full = try support.loadIndex(allocator, io, paths.fasta_path, .on_access);
    defer full.deinit();
    var selected = try support.loadIndex(allocator, io, paths.fasta_path, .{ .requested = &.{ "a", "b", "missing" } });
    defer selected.deinit();
    try std.testing.expectEqual(@as(usize, 2), selected.records.len());
    try std.testing.expectEqualStrings("a", selected.recordName(0).?);
    try std.testing.expectEqualStrings("b", selected.recordName(1).?);
    try std.testing.expectEqual(full.records.len(), selected.records.len());
    for (0..full.records.len()) |i| {
        const actual = try recordAt(&selected, i);
        const expected = try recordAt(&full, i);
        try std.testing.expectEqual(expected.seq_offset, actual.seq_offset);
        try std.testing.expectEqual(expected.seq_len, actual.seq_len);
    }
    try std.testing.checkAllAllocationFailures(allocator, loadAndDeinitForAllocationCheck, .{ paths.fasta_path, main.index_format.LoadMode{ .requested = &.{ "a", "b" } } });
    var scanned = try support.loadIndex(allocator, io, paths.fasta_path, .stats_scan);
    defer scanned.deinit();
    try std.testing.expectEqual(full.records.len(), scanned.records.len());
    var single = try support.loadIndex(allocator, io, paths.fasta_path, .{ .requested = &.{"a"} });
    defer single.deinit();
    try std.testing.expectEqual(@as(usize, 1), single.records.len());
    try std.testing.expectEqual((try recordAt(&full, 0)).seq_len, (try recordAt(&single, 0)).seq_len);
}

test "[regression] - [fai loader]: empty numeric columns never shift later values" {
    const allocator = std.testing.allocator;
    const rows = [_][]const u8{
        "seq1\t\t4\t6\t4\t5",
        "seq1\t4\t\t6\t4\t5",
        "seq1\t4\t6\t\t4\t5",
        "seq1\t4\t6\t4\t\t5",
    };
    const modes = [_]main.index_format.LoadMode{
        .catalog,
        .{ .requested = &.{"seq1"} },
        .stats_scan,
    };
    for (rows) |row| {
        for ([_][]const u8{ "", "\n" }) |ending| {
            const fai = try std.mem.concat(allocator, u8, &.{ row, ending });
            defer allocator.free(fai);
            const paths = try writeFastaAndFai(allocator, "empty-fai-column", ">seq1\nACGT\n", fai);
            defer allocator.free(paths.fasta_path);
            defer allocator.free(paths.fai_path);
            defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
            defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
            for (modes) |mode| {
                try expectLoadFailure(allocator, paths.fasta_path, mode, .@"index.fai.field_missing");
            }
        }
    }
}

test "[regression] - [fai loader]: long rows preserve empty numeric columns" {
    const allocator = std.testing.allocator;
    const rows = [_][]const u8{
        "a\t\t4\t3\t4\t5\t",
        "a\t4\t\t3\t4\t5\t",
        "a\t4\t3\t\t4\t5\t",
        "a\t4\t3\t4\t\t5\t",
    };
    for (rows) |row| {
        for ([_][]const u8{ "", "\n" }) |ending| {
            var text = std.Io.Writer.Allocating.init(allocator);
            defer text.deinit();
            try text.writer.writeAll(row);
            try text.writer.splatByteAll('x', 70_000);
            try text.writer.writeAll(ending);
            const paths = try writeFastaAndFai(allocator, "empty-long-fai-column", ">a\nACGT\n", text.written());
            defer allocator.free(paths.fasta_path);
            defer allocator.free(paths.fai_path);
            defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
            defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
            for ([_]main.index_format.LoadMode{ .catalog, .{ .requested = &.{"a"} }, .{ .requested = &.{"missing"} } }) |mode| {
                try expectLoadFailure(allocator, paths.fasta_path, mode, .@"index.fai.field_missing");
            }
        }
    }
}

test "[regression] - [fai loader]: valid column positions retain decimal and trailing-column acceptance" {
    const allocator = std.testing.allocator;
    const cases = .{
        .{ ">seq1\nACGT\n", "seq1", "seq1\t0004\t0006\t0004\t0005\tignored\t", @as(u64, 4), @as(u64, 6), @as(u32, 4), @as(u32, 5) },
        .{ ">\nACGT\n", "", "\t4\t2\t4\t5\n", @as(u64, 4), @as(u64, 2), @as(u32, 4), @as(u32, 5) },
        .{ ">s\nA\n", "s", "s\t1\t3\t4294967295\t4294967295\n", @as(u64, 1), @as(u64, 3), std.math.maxInt(u32), std.math.maxInt(u32) },
    };
    inline for (cases) |case| {
        const paths = try writeFastaAndFai(allocator, "valid-fai-columns", case[0], case[2]);
        defer allocator.free(paths.fasta_path);
        defer allocator.free(paths.fai_path);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        const names = [_][]const u8{case[1]};
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &names }, .stats_scan }) |mode| {
            var index = try support.loadIndex(allocator, io, paths.fasta_path, mode);
            defer index.deinit();
            try std.testing.expectEqual(@as(usize, 1), index.records.len());
            if (mode != .stats_scan) try std.testing.expectEqualStrings(case[1], index.recordName(0).?);
            try std.testing.expectEqual(case[1].len, (try recordAt(&index, 0)).name_len);
            try std.testing.expectEqual(case[3], (try recordAt(&index, 0)).seq_len);
            try std.testing.expectEqual(case[4], (try recordAt(&index, 0)).seq_offset);
            try std.testing.expectEqual(case[5], (try recordAt(&index, 0)).line_bases);
            try std.testing.expectEqual(case[6], (try recordAt(&index, 0)).line_bytes);
        }
    }
}

test "[regression] - [fai loader]: get and stats reject shifted empty columns" {
    const allocator = std.testing.allocator;
    const paths = try writeFastaAndFai(allocator, "shifted-fai-cli", ">seq1\nACGT\n", "seq1\t\t4\t6\t4\t5\n");
    defer allocator.free(paths.fasta_path);
    defer allocator.free(paths.fai_path);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
    const names = try writeFastaArtifact(allocator, "shifted-fai-names", "seq1\n");
    defer allocator.free(names);
    defer std.Io.Dir.cwd().deleteFile(io, names) catch {};
    const bed = try writeFastaArtifact(allocator, "shifted-fai-bed", "seq1\t0\t4\n");
    defer allocator.free(bed);
    defer std.Io.Dir.cwd().deleteFile(io, bed) catch {};
    try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", paths.fasta_path, "seq1" }, .@"index.fai.field_missing");
    try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", paths.fasta_path, "--names", names }, .@"index.fai.field_missing");
    try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "get", paths.fasta_path, "--bed", bed }, .@"index.fai.field_missing");
    try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "stats", paths.fasta_path }, .@"index.fai.field_missing");
}

test "[failure] - [index diagnostics]: FAI fields retain original rows and byte columns" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const Case = struct { row: []const u8, code: diagnostics.Code, field: diagnostics.FaiField, column: u64, value: ?[]const u8 = null };
    const cases = [_]Case{
        .{ .row = "seq1", .code = .@"index.fai.field_missing", .field = .length, .column = 5 },
        .{ .row = "seq1\t\t4\t6\t4\t5", .code = .@"index.fai.field_missing", .field = .length, .column = 6 },
        .{ .row = "seq1\t4\t\t6\t4\t5", .code = .@"index.fai.field_missing", .field = .offset, .column = 8 },
        .{ .row = "seq1\t4\t6\t\t4\t5", .code = .@"index.fai.field_missing", .field = .line_bases, .column = 10 },
        .{ .row = "seq1\t4\t6\t4\t\t5", .code = .@"index.fai.field_missing", .field = .line_bytes, .column = 12 },
        .{ .row = "seq1\t4\t6\t4", .code = .@"index.fai.field_missing", .field = .line_bytes, .column = 11 },
        .{ .row = "seq1\t18446744073709551616\t6\t4\t5", .code = .@"index.fai.integer_out_of_range", .field = .length, .column = 6, .value = "18446744073709551616" },
        .{ .row = "seq1\t4\t18446744073709551616\t4\t5", .code = .@"index.fai.integer_out_of_range", .field = .offset, .column = 8, .value = "18446744073709551616" },
        .{ .row = "seq1\t4\t6\t4294967296\t5", .code = .@"index.fai.integer_out_of_range", .field = .line_bases, .column = 10, .value = "4294967296" },
        .{ .row = "seq1\t4\t6\t4\t4294967296", .code = .@"index.fai.integer_out_of_range", .field = .line_bytes, .column = 12, .value = "4294967296" },
        .{ .row = "seq1\t18446744073709551616x\t6\t4\t5", .code = .@"index.fai.invalid_integer", .field = .length, .column = 6, .value = "18446744073709551616x" },
        .{ .row = "seq1\t4\t-6\t4\t5", .code = .@"index.fai.invalid_integer", .field = .offset, .column = 8, .value = "-6" },
    };
    const prefix = "\nseq1\t4\t6\t4\t5\n\n";
    for (cases) |case| {
        const text = try std.mem.concat(allocator, u8, &.{ prefix, case.row });
        const paths = try writeFastaAndFai(allocator, "fai-diagnostic-fields", ">seq1\nACGT\n", text);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        for ([_]main.index_format.LoadMode{ .catalog, .stats_scan }) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, case.code, &failure);
            const primary = failure.primary().?;
            try std.testing.expectEqual(diagnostics.SourceId.fai, primary.at.source.?);
            try std.testing.expectEqual(@as(u64, 4), primary.at.location.line.?);
            try std.testing.expectEqual(@as(u64, 2), primary.at.record.?);
            try std.testing.expectEqual(case.column, primary.at.location.column.?);
            try std.testing.expectEqual(prefix.len + case.column - 1, primary.at.location.offset.?);
            switch (primary.details) {
                .@"index.fai.field_missing" => |facts| {
                    try std.testing.expectEqual(case.field, facts.field);
                    try std.testing.expectEqual(@intFromEnum(case.field) + 1, facts.column_number);
                },
                .@"index.fai.invalid_integer" => |facts| {
                    try std.testing.expectEqual(case.field, facts.field);
                    try std.testing.expectEqualStrings(case.value.?, failure.preview(facts.value.?));
                },
                .@"index.fai.integer_out_of_range" => |facts| {
                    try std.testing.expectEqual(case.field, facts.field);
                    try std.testing.expectEqualStrings(case.value.?, failure.preview(facts.value.?));
                    try std.testing.expectEqual(if (case.field == .length or case.field == .offset) @as(u64, std.math.maxInt(u64)) else @as(u64, std.math.maxInt(u32)), facts.max_value);
                },
                else => unreachable,
            }
            try std.testing.expectEqualStrings(paths.fasta_path, failure.sourceLabel(.fai).?);
            try std.testing.expectEqual(diagnostics.Suffix.fai, failure.sources.get(.fai).?.suffix);
        }
    }
}

test "[failure] - [index diagnostics]: FAI sequence extents retain wide products and remainders" {
    const allocator = std.testing.allocator;
    const cases = [_]struct { row: []const u8, required_end: u128 }{
        .{ .row = "a\t18446744073709551615\t3\t1\t4294967295\n", .required_end = 79228162495817593511244464134 },
        .{ .row = "a\t18446744073709551615\t3\t4\t5\n", .required_end = 23058430092136939521 },
        .{ .row = "a\t8\t3\t4\t5\n", .required_end = 12 },
    };
    for (cases) |case| {
        const paths = try writeFastaAndFai(allocator, "fai-wide-extent", ">a\nACGT\n", case.row);
        defer allocator.free(paths.fasta_path);
        defer allocator.free(paths.fai_path);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        for ([_]main.index_format.LoadMode{ .catalog, .{ .requested = &.{"a"} }, .stats_scan }) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.sequence_span_out_of_bounds", &failure);
            const facts = failure.primary().?.details.@"index.sequence_span_out_of_bounds";
            try std.testing.expectEqual(case.required_end, facts.required_end);
            try std.testing.expectEqual(@as(u64, 8), facts.source_bytes);
            try std.testing.expectEqual(@as(u64, 1), failure.primary().?.at.record.?);
        }
    }
}

test "[failure] - [index diagnostics]: long FAI fields retain original lengths and lexical priority" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_]bool{ false, true }) |invalid_byte| {
        var text = std.Io.Writer.Allocating.init(allocator);
        defer text.deinit();
        try text.writer.writeAll("a\t4\t3\t");
        try text.writer.splatByteAll('9', 70_000);
        if (invalid_byte) try text.writer.writeByte('x');
        try text.writer.writeAll("\t5\textra\n");
        const paths = try writeFastaAndFai(allocator, "long-fai-diagnostic", ">a\nACGT\n", text.written());
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.fai_path) catch {};
        const code: diagnostics.Code = if (invalid_byte) .@"index.fai.invalid_integer" else .@"index.fai.integer_out_of_range";
        for ([_]main.index_format.LoadMode{ .catalog, .stats_scan, .{ .requested = &.{"missing"} }, .{ .requested = &.{"a"} } }) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, code, &failure);
            const primary = failure.primary().?;
            const value = if (invalid_byte) primary.details.@"index.fai.invalid_integer".value.? else primary.details.@"index.fai.integer_out_of_range".value.?;
            try std.testing.expectEqual(@as(u64, 70_000) + @intFromBool(invalid_byte), value.original_bytes);
            try std.testing.expectEqual(@as(u16, 32), value.len);
            try std.testing.expect(value.truncated);
            try std.testing.expectEqual(@as(u64, 6), primary.at.location.offset.?);
            try std.testing.expectEqual(@as(u64, 7), primary.at.location.column.?);
        }
    }
}

test "[failure] - [index diagnostics]: ZFI header and record facts agree in all modes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const original = try zfiBytes(allocator, &index, fasta.len, 0);
    const Mutation = struct { offset: usize, width: u8, value: u64, code: diagnostics.Code };
    const cases = [_]Mutation{
        .{ .offset = 0, .width = 1, .value = 'X', .code = .@"index.invalid_magic" },
        .{ .offset = 3, .width = 1, .value = 255, .code = .@"index.unsupported_version" },
        .{ .offset = 4, .width = 4, .value = 0, .code = .@"index.no_records" },
        .{ .offset = 4, .width = 4, .value = std.math.maxInt(u32), .code = .@"index.truncated" },
        .{ .offset = 8, .width = 8, .value = 999, .code = .@"index.source_size_mismatch" },
        .{ .offset = 30, .width = 2, .value = 0xFFFF, .code = .@"index.name_out_of_bounds" },
        .{ .offset = 22, .width = 2, .value = 2, .code = .@"index.name_out_of_bounds" },
        .{ .offset = 16, .width = 4, .value = fasta.len, .code = .@"index.sequence_offset_out_of_bounds" },
        .{ .offset = 24, .width = 4, .value = 0, .code = .@"index.zero_sequence_length" },
        .{ .offset = 24, .width = 8, .value = main.index_format.MAX_ZFI_SOURCE_OFFSET, .code = .@"index.sequence_span_out_of_bounds" },
        .{ .offset = 35, .width = 4, .value = 0, .code = .@"index.invalid_line_geometry" },
        .{ .offset = original.len - 12, .width = 1, .value = 'X', .code = .@"index.name_footer_missing" },
        .{ .offset = original.len - 24, .width = 1, .value = 'X', .code = .@"index.name_footer_missing" },
        .{ .offset = original.len - 8, .width = 8, .value = std.math.maxInt(u64), .code = .@"index.name_blob_out_of_bounds" },
    };
    for (cases) |case| {
        const bytes = try allocator.dupe(u8, original);
        switch (case.width) {
            1 => bytes[case.offset] = @intCast(case.value),
            2 => std.mem.writeInt(u16, bytes[case.offset..][0..2], @intCast(case.value), .little),
            4 => std.mem.writeInt(u32, bytes[case.offset..][0..4], @intCast(case.value), .little),
            8 => std.mem.writeInt(u64, bytes[case.offset..][0..8], case.value, .little),
            else => unreachable,
        }
        const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-record", fasta, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        var first: diagnostics.Diagnostic = undefined;
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"a"} }, .stats_scan, .{ .requested = &.{"missing"} } }, 0..) |mode, i| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, case.code, &failure);
            const primary = failure.primary().?;
            try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source.?);
            if (i == 0) first = primary.* else {
                try std.testing.expectEqualDeep(first.at, primary.at);
                try std.testing.expectEqualDeep(first.details, primary.details);
            }
            try std.testing.expect(primary.cause == null);
            switch (primary.details) {
                .@"index.invalid_magic" => |facts| {
                    try std.testing.expectEqual(ZFI_MAGIC, facts.expected_magic);
                    try std.testing.expectEqualStrings("XFI\x02", failure.preview(facts.actual.?));
                },
                .@"index.unsupported_version" => |facts| {
                    try std.testing.expectEqual(@as(u8, 255), facts.found_version);
                    try std.testing.expectEqual(@as(u8, 2), facts.supported_version);
                },
                .@"index.no_records" => |facts| try std.testing.expectEqual(diagnostics.Format.zfi, facts.format),
                .@"index.truncated" => |facts| {
                    try std.testing.expectEqual(diagnostics.Section.records, facts.section);
                    try std.testing.expectEqual(@as(u128, 16) + @as(u128, std.math.maxInt(u32)) * 24, facts.required_end);
                    try std.testing.expectEqual(bytes.len, facts.available_end);
                },
                .@"index.source_size_mismatch" => |facts| {
                    try std.testing.expectEqual(@as(u64, 999), facts.stored_bytes);
                    try std.testing.expectEqual(fasta.len, facts.actual_bytes);
                },
                .@"index.name_out_of_bounds" => |facts| {
                    try std.testing.expectEqual(@as(u64, 1), facts.blob_bytes);
                    try std.testing.expectEqual(if (case.offset == 30) @as(u64, 0xFFFF) else @as(u64, 0), facts.name_offset);
                    try std.testing.expectEqual(if (case.offset == 22) @as(u64, 2) else @as(u64, 1), facts.name_bytes);
                },
                .@"index.zero_sequence_length" => try std.testing.expectEqual(@as(u64, 1), primary.at.record.?),
                .@"index.sequence_offset_out_of_bounds" => |facts| {
                    try std.testing.expectEqual(fasta.len, facts.sequence_offset);
                    try std.testing.expectEqual(fasta.len, facts.source_bytes);
                },
                .@"index.sequence_span_out_of_bounds" => |facts| {
                    try std.testing.expectEqual(@as(u128, 5 << 46) + 1, facts.required_end);
                    try std.testing.expectEqual(fasta.len, facts.source_bytes);
                    try std.testing.expect(facts.entry == null);
                },
                .@"index.invalid_line_geometry" => |facts| {
                    try std.testing.expectEqual(@as(u64, 0), facts.geometry.line_bases);
                    try std.testing.expectEqual(@as(u64, 1), facts.geometry.line_bytes);
                    try std.testing.expectEqual(diagnostics.GeometryConstraint.bases_zero, facts.constraint);
                },
                .@"index.name_footer_missing" => {},
                .@"index.name_blob_out_of_bounds" => |facts| {
                    try std.testing.expectEqual(std.math.maxInt(u64), facts.name_blob_bytes);
                    try std.testing.expectEqual(@as(u64, 41), facts.payload_end);
                    try std.testing.expectEqual(@as(u64, 40), facts.records_end);
                },
                else => return error.UnexpectedDiagnostic,
            }
        }
    }
}

test "[failure] - [index diagnostics]: side tables retain wire offsets and operands" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nAAA\nCCCC\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const original = try zfiBytes(allocator, &index, fasta.len, 0);
    const side: usize = @intCast(index.records.items[0].sideTableOffset());
    const first = side + 8;
    const second = first + main.index_format.ZFI_SIDE_TABLE_LINE_BYTES;
    const Case = struct { offset: usize, width: u8 = 8, value: u64, code: diagnostics.Code };
    const cases = [_]Case{
        .{ .offset = side, .value = 0, .code = .@"index.side_table_empty" },
        .{ .offset = side, .value = std.math.maxInt(u64), .code = .@"index.arithmetic_overflow" },
        .{ .offset = side, .value = 100, .code = .@"index.truncated" },
        .{ .offset = first, .value = 1, .code = .@"index.side_base_start_mismatch" },
        .{ .offset = first + 20, .width = 4, .value = 0, .code = .@"index.side_line_geometry_invalid" },
        .{ .offset = first + 16, .width = 4, .value = 0, .code = .@"index.side_line_geometry_invalid" },
        .{ .offset = first + 16, .width = 4, .value = 2, .code = .@"index.side_line_geometry_invalid" },
        .{ .offset = first + 8, .value = std.math.maxInt(u64), .code = .@"index.arithmetic_overflow" },
        .{ .offset = first + 16, .width = 4, .value = fasta.len, .code = .@"index.sequence_span_out_of_bounds" },
        .{ .offset = first + 8, .value = 4, .code = .@"index.side_first_offset_mismatch" },
        .{ .offset = second + 8, .value = 3, .code = .@"index.side_line_order_invalid" },
        .{ .offset = 24, .value = 10, .code = .@"index.side_sequence_length_mismatch" },
    };
    for (cases) |case| {
        const bytes = try allocator.dupe(u8, original);
        if (case.width == 4) std.mem.writeInt(u32, bytes[case.offset..][0..4], @intCast(case.value), .little) else std.mem.writeInt(u64, bytes[case.offset..][0..8], case.value, .little);
        const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-side", fasta, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        var expected: diagnostics.Diagnostic = undefined;
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"a"} }, .{ .requested = &.{"missing"} } }, 0..) |mode, i| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, case.code, &failure);
            const primary = failure.primary().?;
            try std.testing.expectEqual(@as(u64, 1), primary.at.record.?);
            try std.testing.expect(primary.at.location.offset.? >= side);
            if (i == 0) expected = primary.* else {
                try std.testing.expectEqualDeep(expected.at, primary.at);
                try std.testing.expectEqualDeep(expected.details, primary.details);
            }
            try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source.?);
            try std.testing.expectEqual(diagnostics.Stage.parse, primary.stage);
            try std.testing.expect(primary.cause == null);
            switch (primary.details) {
                .@"index.side_table_empty" => try std.testing.expectEqual(side, primary.at.location.offset.?),
                .@"index.arithmetic_overflow" => |facts| {
                    try std.testing.expectEqual(case.value, facts.lhs);
                    try std.testing.expectEqual(if (case.offset == side) @as(u128, 24) else @as(u128, 4), facts.rhs);
                    try std.testing.expectEqual(if (case.offset == side) diagnostics.IndexField.side_count else diagnostics.IndexField.side_byte_end, facts.field);
                    try std.testing.expectEqual(if (case.offset == side) diagnostics.Arithmetic.multiply else diagnostics.Arithmetic.add, facts.arithmetic);
                    try std.testing.expectEqual(std.math.maxInt(usize), facts.max_value);
                },
                .@"index.truncated" => |facts| {
                    try std.testing.expectEqual(diagnostics.Section.side_table_entries, facts.section);
                    try std.testing.expectEqual(side + 8 + 100 * 24, facts.required_end);
                    try std.testing.expectEqual(side + 8 + 3 * 24, facts.available_end);
                },
                .@"index.side_base_start_mismatch" => |facts| {
                    try std.testing.expectEqual(@as(u64, 1), facts.entry);
                    try std.testing.expectEqual(@as(u64, 1), facts.base_start);
                    try std.testing.expectEqual(@as(u64, 0), facts.expected_base_start);
                },
                .@"index.side_line_geometry_invalid" => |facts| {
                    try std.testing.expectEqual(@as(u64, 1), facts.entry);
                    try std.testing.expectEqual(if (case.offset == first + 20) case.value else @as(u64, 3), facts.geometry.line_bases);
                    try std.testing.expectEqual(if (case.offset == first + 16) case.value else @as(u64, 4), facts.geometry.line_bytes);
                    try std.testing.expectEqual(if (case.offset == first + 20) diagnostics.GeometryConstraint.bases_zero else if (case.value == 0) diagnostics.GeometryConstraint.bytes_zero else diagnostics.GeometryConstraint.bytes_below_bases, facts.constraint);
                },
                .@"index.sequence_span_out_of_bounds" => |facts| {
                    try std.testing.expectEqual(fasta.len + 3, facts.required_end);
                    try std.testing.expectEqual(fasta.len, facts.source_bytes);
                    try std.testing.expectEqual(@as(?u64, 1), facts.entry);
                },
                .@"index.side_first_offset_mismatch" => |facts| {
                    try std.testing.expectEqual(@as(u64, 4), facts.first_byte_offset);
                    try std.testing.expectEqual(@as(u64, 3), facts.sequence_offset);
                },
                .@"index.side_line_order_invalid" => |facts| {
                    try std.testing.expectEqual(@as(u64, 2), facts.entry);
                    try std.testing.expectEqual(@as(u64, 3), facts.byte_offset);
                    try std.testing.expectEqual(@as(u64, 7), facts.minimum_offset);
                },
                .@"index.side_sequence_length_mismatch" => |facts| {
                    try std.testing.expectEqual(@as(u64, 10), facts.indexed_length);
                    try std.testing.expectEqual(@as(u64, 9), facts.side_table_bases);
                },
                else => return error.UnexpectedDiagnostic,
            }
        }
    }
    for ([_]u64{ 16, side + 1 }) |offset| {
        var rec = index.records.items[0];
        try rec.markNonUniform(offset);
        const bytes = try allocator.dupe(u8, original);
        @memcpy(bytes[16..][0..24], std.mem.asBytes(&main.index_format.ZfiRecord.encode(rec)));
        const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-side-offset", fasta, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"missing"} } }) |mode| {
            var failure: diagnostics.Failure = undefined;
            const code: diagnostics.Code = if (offset == 16) .@"index.side_table_offset_out_of_bounds" else .@"index.side_table_misaligned";
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, code, &failure);
            if (offset == 16) try std.testing.expectEqual(offset, failure.primary().?.details.@"index.side_table_offset_out_of_bounds".offset) else try std.testing.expectEqual(offset, failure.primary().?.details.@"index.side_table_misaligned".offset);
        }
    }
}

test "[regression] - [index diagnostics]: trailer defects precede record checks" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n>b\nTTTT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const original = try zfiBytes(allocator, &index, fasta.len, 0);
    for ([_]enum { footer, identity, blob }{ .footer, .identity, .blob }) |defect| {
        const bytes = try allocator.dupe(u8, original);
        std.mem.writeInt(u32, bytes[24..][0..4], 0, .little);
        switch (defect) {
            .footer => bytes[bytes.len - 12] = 'X',
            .identity => bytes[bytes.len - 24] = 'X',
            .blob => std.mem.writeInt(u64, bytes[bytes.len - 8 ..][0..8], std.math.maxInt(u64), .little),
        }
        const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-priority", fasta, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        const code: diagnostics.Code = if (defect == .blob) .@"index.name_blob_out_of_bounds" else .@"index.name_footer_missing";
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"missing"} } }) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, code, &failure);
            const primary = failure.primary().?;
            try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source.?);
            try std.testing.expectEqual(diagnostics.Stage.parse, primary.stage);
            try std.testing.expectEqual(bytes.len, primary.at.location.offset.?);
            try std.testing.expect(primary.at.record == null);
            try std.testing.expectEqual(@as(usize, 0), failure.secondaryItems().len);
            if (defect == .blob) {
                const facts = primary.details.@"index.name_blob_out_of_bounds";
                try std.testing.expectEqual(bytes.len - main.index_format.ZFI_NAME_FOOTER_BYTES - main.index_format.ZFI_SOURCE_ID_BYTES - 4 * main.index_format.NameTable.slotCount(2), facts.payload_end);
                try std.testing.expectEqual(main.index_format.ZFI_HEADER_BYTES + 2 * main.index_format.ZFI_INDEX_RECORD_BYTES, facts.records_end);
            }
        }
    }
}

test "[failure] - [index diagnostics]: required reads and stats retain causes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-read-diagnostic", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    for ([_]bool{ false, true }) |short_read| {
        var vtable = io.vtable.*;
        vtable.fileReadPositional = if (short_read) shortIndexRead else failIndexRead;
        const faulty_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
        var failure: diagnostics.Failure = undefined;
        failure.reset(.get);
        var stream: main.index_format.EntryStream = undefined;
        try std.testing.expectError(error.Io, stream.open(allocator, faulty_io, paths.fasta_path, .automatic, &failure));
        const primary = failure.primary().?;
        try std.testing.expectEqual(diagnostics.SourceId.zfi, primary.at.source.?);
        try std.testing.expectEqual(diagnostics.Stage.read, primary.stage);
        if (short_read) {
            const facts = primary.details.@"io.unexpected_eof".exact;
            try std.testing.expectEqual(@as(u64, 16), facts.requested_bytes);
            try std.testing.expectEqual(@as(u64, 0), facts.received_bytes);
            try std.testing.expectEqual(@as(u64, 0), primary.at.location.offset.?);
        } else {
            try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", primary.code());
            try std.testing.expectEqual(@as(anyerror, error.InputOutput), primary.cause.?.zig);
        }
    }
    var vtable = io.vtable.*;
    vtable.fileStat = failIndexStat;
    const faulty_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectError(error.Io, main.index_format.loadIndex(allocator, faulty_io, paths.fasta_path, .automatic, .lookup_full_map, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.stat_failed", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.SourceId.fasta, failure.primary().?.at.source.?);
    try std.testing.expectEqual(@as(anyerror, error.AccessDenied), failure.primary().?.cause.?.zig);
    failure.reset(.stats);
    var loaded = try main.index_format.loadIndex(allocator, io, paths.fasta_path, .automatic, .on_access, &failure);
    defer loaded.deinit();
    try std.testing.expectEqual(diagnostics.Operation.stats, failure.operation);
    try std.testing.expectEqual(null, failure.primary());
}

test "[failure] - [index diagnostics]: discovery and allocation failures remain recoverable" {
    const allocator = std.testing.allocator;
    const fasta_path = try writeFastaArtifact(allocator, "index-diagnostic-discovery", ">a\nACGT\n");
    defer allocator.free(fasta_path);
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    var failure: diagnostics.Failure = undefined;
    try expectIndexLoadCode(allocator, fasta_path, .lookup_full_map, .@"index.not_found", &failure);
    try std.testing.expectEqual(diagnostics.Stage.open, failure.primary().?.stage);
    try std.testing.expectEqual(@as(usize, 0), failure.secondaryItems().len);
    try std.testing.expectEqual(diagnostics.SourceId.zfi, failure.primary().?.details.@"index.not_found".zfi_source);
    const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{fasta_path});
    defer allocator.free(fai_path);
    defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};
    const file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
    try std.Io.File.writeStreamingAll(file, io, "a\t4\t3\t4\t5\n");
    file.close(io);
    for ([_]main.index_format.LoadMode{ .lookup_full_map, .stats_scan, .{ .requested = &.{"a"} } }) |mode| {
        var failing = std.testing.FailingAllocator.init(allocator, .{ .fail_index = 0 });
        try expectIndexLoadCode(failing.allocator(), fasta_path, mode, .@"resource.out_of_memory", &failure);
        try std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage);
        try std.testing.expectEqual(@as(anyerror, error.OutOfMemory), failure.primary().?.cause.?.zig);
        try std.testing.expectEqual(failing.allocated_bytes, failing.freed_bytes);
    }
    try expectIndexLoadCode(allocator, "zig-cache/test-artifacts/missing-index-diagnostic.fa", .lookup_full_map, .@"io.file_not_found", &failure);
    try std.testing.expectEqual(diagnostics.SourceId.fasta, failure.primary().?.at.source.?);
    try std.testing.expectEqual(diagnostics.Stage.open, failure.primary().?.stage);
}

test "[failure] - [index diagnostics]: truncation and source identity precede record geometry" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    defer index.deinit(allocator);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    for ([_]usize{ 0, 4, 15, 20 }) |len| {
        const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-truncation", fasta, bytes[0..len]);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"a"} } }) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.truncated", &failure);
            const facts = failure.primary().?.details.@"index.truncated";
            try std.testing.expectEqual(len, facts.available_end);
            try std.testing.expectEqual(if (len < 16) @as(u128, 16) else @as(u128, 40), facts.required_end);
        }
    }
    const paths = try writeFastaAndRawZfi(allocator, "zfi-diagnostic-identity", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    {
        const file = try std.Io.Dir.cwd().openFile(io, paths.zfi_path, .{ .mode = .read_write });
        defer file.close(io);
        const zeros: [8]u8 = .{0} ** 8;
        try file.writePositionalAll(io, &zeros, 24);
        try file.writePositionalAll(io, &zeros, bytes.len - 20);
    }
    for ([_]main.index_format.LoadMode{ .lookup_full_map, .{ .requested = &.{"a"} } }) |mode| {
        var failure: diagnostics.Failure = undefined;
        try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.source_mtime_mismatch", &failure);
        const primary = failure.primary().?;
        try std.testing.expectEqual(@as(u64, 0), primary.details.@"index.source_mtime_mismatch".stored_mtime_ns);
        try std.testing.expectEqual(diagnostics.Stage.check_identity, primary.stage);
        try std.testing.expectEqual(bytes.len - 20, primary.at.location.offset.?);
        try std.testing.expectEqual(@as(i128, try statMtimeNs(paths.fasta_path)), primary.details.@"index.source_mtime_mismatch".actual_mtime_ns);
    }
}

test "[failure] - [index diagnostics]: retains source, write, publish and comparison causes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "index-diagnostics", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    const Case = struct { code: diagnostics.Code, stage: diagnostics.Stage, cause: anyerror };
    const cases = [_]Case{
        .{ .code = .@"io.stat_failed", .stage = .stat, .cause = error.AccessDenied },
        .{ .code = .@"io.read_failed", .stage = .read, .cause = error.InputOutput },
        .{ .code = .@"io.create_failed", .stage = .create, .cause = error.AccessDenied },
        .{ .code = .@"io.flush_failed", .stage = .flush, .cause = error.NoSpaceLeft },
        .{ .code = .@"io.publish_failed", .stage = .publish, .cause = error.AccessDenied },
        .{ .code = .@"input.changed", .stage = .compare_source, .cause = error.FileNotFound },
    };
    for ([_]bool{ false, true }) |large| {
        if (large) {
            const source = try std.Io.Dir.cwd().openFile(io, path, .{ .mode = .read_write });
            defer source.close(io);
            const wraps = "ACGT\n" ** 2048;
            var offset: u64 = 10;
            while (offset <= 32 * 1024 * 1024) : (offset += wraps.len) try source.writePositionalAll(io, wraps, offset);
        }
        const source_bytes = (try std.Io.Dir.cwd().statFile(io, path, .{})).size;
        for (cases, 0..) |case, i| {
            try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = zfi_path, .data = "previous index" });
            var vtable = io.vtable.*;
            switch (i) {
                0 => vtable.fileStat = failIndexStat,
                1 => vtable.fileReadPositional = failIndexRead,
                2 => vtable.dirCreateFileAtomic = failIndexCreate,
                3 => vtable.fileWritePositional = failIndexWrite,
                4 => vtable.dirRename = failIndexPublish,
                5 => vtable.dirStatFile = disappearIndexSource,
                else => unreachable,
            }
            var failure: diagnostics.Failure = undefined;
            try expectIndexRunCode(.{ .userdata = io.userdata, .vtable = &vtable }, path, false, if (large and i == 3) .@"io.write_failed" else case.code, &failure);
            const primary = failure.primary().?;
            try std.testing.expectEqual(if (large and i == 3) diagnostics.Stage.write else case.stage, primary.stage);
            try std.testing.expectEqual(case.cause, primary.cause.?.zig);
            try std.testing.expectEqualStrings("previous index", try readTestFile(allocator, zfi_path));
            try std.testing.expectEqualStrings(path, failure.sourceLabel(.zfi).?);
            if (i == 3 or i == 5) try std.testing.expectEqual(.aborted, failure.outputs.file.state);
            if (i == 4) try std.testing.expectEqual(.unknown, failure.outputs.file.state);
            if (i == 5) {
                const changed = primary.details.@"input.changed";
                try std.testing.expect(changed.disappeared);
                try std.testing.expectEqual(@as(?u64, source_bytes), changed.before.bytes);
                try std.testing.expectEqual(@as(?diagnostics.FileSnapshot, null), changed.after);
            }
        }
    }
}

test "[failure] - [index]: every allocation failure is reported and released" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_][]const u8{ ">a\nACGT\nAC\n>b\nGG\nG\nGG\n", ">a\nACGT\nAC\n>b\nGGG\n" }, [_]bool{ false, true }) |data, emit_fai| {
        const path = try writeFastaArtifact(allocator, "index-allocation", data);
        defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
        const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
        defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
        try std.testing.checkAllAllocationFailures(std.testing.allocator, runIndexForAllocationCheck, .{ path, emit_fai });
    }
}

test "[failure] - [index]: a large input reports and releases its first allocation failures" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const line = "ACGT" ** 15 ++ "\n";
    var data: std.ArrayList(u8) = .empty;
    try data.appendSlice(allocator, ">large\n");
    while (data.items.len <= 33 * 1024 * 1024) try data.appendSlice(allocator, line);
    const path = try writeFastaArtifact(allocator, "index-large-allocation", data.items);
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    for (0..4) |fail_index| {
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = fail_index });
        try std.testing.expectError(error.OutOfMemory, runIndexForAllocationCheck(failing.allocator(), path, false));
        try std.testing.expect(failing.has_induced_failure);
    }
}

test "[failure] - [index diagnostics]: stages complete output and retains cleanup secondaries" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "fai-diagnostics", ">seq\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    for (0..5) |i| {
        var vtable = io.vtable.*;
        vtable.operate = discardIndexStdout;
        switch (i) {
            0 => vtable.dirCreateFile = failIndexSpoolCreate,
            1 => vtable.fileWritePositional = failIndexWrite,
            2 => vtable.operate = failIndexStdout,
            3 => vtable.dirDeleteFile = failIndexCleanup,
            4 => {
                vtable.fileWritePositional = failIndexWrite;
                vtable.dirDeleteFile = failIndexCleanup;
            },
            else => unreachable,
        }
        const code: diagnostics.Code = switch (i) {
            0 => .@"io.temp_unavailable",
            3 => .@"io.cleanup_failed",
            else => .@"io.flush_failed",
        };
        var failure: diagnostics.Failure = undefined;
        try expectIndexRunCode(.{ .userdata = io.userdata, .vtable = &vtable }, path, true, code, &failure);
        const primary = failure.primary().?;
        if (i == 0) {
            try std.testing.expectEqual(diagnostics.TemporaryPurpose.fai, primary.details.@"io.temp_unavailable".purpose);
            try std.testing.expect(primary.details.@"io.temp_unavailable".attempts.len > 0);
        } else if (i == 2) {
            try std.testing.expectEqual(diagnostics.SourceId.stdout, primary.at.source.?);
            try std.testing.expectEqual(error.BrokenPipe, primary.cause.?.zig);
            try std.testing.expect(failure.outputs.stdout.may_be_partial);
            try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
        } else {
            try std.testing.expectEqual(diagnostics.SourceId.fai_spool, primary.at.source.?);
            if (i == 3) {
                try std.testing.expect(failure.work_complete);
                try std.testing.expectEqual(.succeeded, failure.outputs.stdout.final_flush);
                try std.testing.expectEqual(@as(?u64, 1), failure.outputs.stdout.records);
            } else {
                try std.testing.expectEqual(error.NoSpaceLeft, primary.cause.?.zig);
                try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
            }
            if (i == 3 or i == 4) {
                defer std.Io.Dir.cwd().deleteFile(io, failure.sourceLabel(.fai_spool).?) catch {};
                if (i == 4) {
                    try std.testing.expectEqual(@as(usize, 1), failure.secondaryItems().len);
                    try std.testing.expectEqual(diagnostics.Code.@"io.cleanup_failed", failure.secondaryItems()[0].code());
                    try std.testing.expectEqual(error.AccessDenied, failure.secondaryItems()[0].cause.?.zig);
                }
            }
        }
    }
}

test "[failure] - [zfi writer]: allocation writer failure is typed OOM" {
    var index = try main.indexer.scanZfiData(std.testing.allocator, ">seq\nACGT\n");
    defer index.deinit(std.testing.allocator);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    try std.testing.expectError(error.OutOfMemory, main.indexer.zfiIndexToBytes(failing.allocator(), &index, 10, 0, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.AllocationRole.serialized_index, failure.primary().?.details.@"resource.out_of_memory".allocation_role);
}

test "[failure] - [index diagnostics]: retains original record for nonuniform FAI rejection" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "nonuniform-diagnostic", ">empty\n>ragged\nACGT\nAC\nTT\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    var failure: diagnostics.Failure = undefined;
    try expectIndexRunCode(io, path, true, .@"index.nonuniform_fai", &failure);
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(?u64, 2), primary.at.record);
    try std.testing.expectEqual(@as(?u64, 7), primary.at.location.offset);
    try std.testing.expectEqual(diagnostics.NonuniformReason.interior_short_line, primary.details.@"index.nonuniform_fai".reason);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
}

test "[integration] - [index diagnostics]: completes ZFI publication before returning success" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "index-result", ">empty\n>seq\nACGT\n>seq\nTT\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    var failure: diagnostics.Failure = undefined;
    const result = try main.indexer.runIndex(std.testing.allocator, io, .empty, path, .{}, &failure);
    try std.testing.expectEqual(@as(u32, 2), result.record_count);
    try std.testing.expectEqual(@as(?*const diagnostics.Diagnostic, null), failure.primary());
    try std.testing.expect(failure.complete and failure.work_complete);
    try std.testing.expectEqual(.published, failure.outputs.file.state);
    try std.testing.expectEqual(.succeeded, failure.outputs.file.final_flush);
    failure.reset(.get);
    var loaded = try main.index_format.loadIndex(allocator, io, path, .automatic, .on_access, &failure);
    defer loaded.deinit();
    try std.testing.expectEqual(@as(usize, 2), loaded.records.len());
    try std.testing.expectEqualStrings("seq", loaded.recordName(0).?);
    try std.testing.expectEqual(@as(?usize, 0), loaded.lookupName("seq"));
}

test "[regression] - [fasta scanner]: empty buffers reject before reading or allocating" {
    const Context = struct {
        fn emit(_: *@This(), _: main.indexer.FastaRecordEmit) !void {
            return error.UnexpectedEmit;
        }
    };
    var context = Context{};
    var buffer: [0]u8 = .{};
    var failure: diagnostics.Failure = undefined;
    failure.reset(.index);
    var reader = std.Io.Reader.fixed(">name\nACGT\n");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });

    try std.testing.expectError(error.InvalidArgument, main.indexer.scanFastaReader(failing.allocator(), &reader, &buffer, .{ .failure = &failure }, &context, Context.emit));
    try std.testing.expectEqual(diagnostics.Code.@"usage.invalid_option_value", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.Argument.read_buffer, failure.primary().?.details.@"usage.invalid_option_value".option);
    try std.testing.expectEqual(diagnostics.Expected.nonempty_buffer, failure.primary().?.details.@"usage.invalid_option_value".expected);
    try std.testing.expectEqual(@as(usize, 0), reader.seek);
    try std.testing.expectEqual(@as(usize, 0), failing.allocations);
    try std.testing.expectError(error.InvalidArgument, main.indexer.scanZfiReader(failing.allocator(), &reader, &buffer, .{ .failure = &failure }));
}

test "[edge] - [index diagnostics]: rejects a FASTA modified before 1970" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta_path = try writeFastaArtifact(allocator, "index-pre-epoch", ">a\nAC\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};
    const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fasta_path});
    defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
    try std.Io.Dir.cwd().setTimestamps(io, fasta_path, .{ .modify_timestamp = .{ .new = .fromNanoseconds(-std.time.ns_per_s) } });

    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.InvalidIndex, main.indexer.runIndex(std.testing.allocator, io, .empty, fasta_path, .{}, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"index.source_timestamp_out_of_range", primary.code());
    const details = primary.details.@"index.source_timestamp_out_of_range";
    try std.testing.expectEqual(@as(i128, -std.time.ns_per_s), details.timestamp_ns);
    try std.testing.expectEqual(@as(i128, 0), details.min_ns);
    try std.testing.expectEqual(@as(i128, std.math.maxInt(u64)), details.max_ns);
    try std.testing.expectError(error.FileNotFound, std.Io.Dir.cwd().statFile(io, zfi_path, .{}));
}

test "[failure] - [index loader]: a FIFO in place of a sidecar fails without waiting for a writer" {
    if (comptime builtin.os.tag != .linux) return error.SkipZigTest;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta_path = try writeFastaArtifact(allocator, "fifo-sidecar", ">a\nACGT\n");
    defer std.Io.Dir.cwd().deleteFile(io, fasta_path) catch {};

    const Case = struct { suffix: []const u8, code: diagnostics.Code };
    for ([_]Case{
        .{ .suffix = "zfi", .code = .@"index.truncated" },
        .{ .suffix = "fai", .code = .@"index.no_records" },
    }) |case| {
        const sidecar = try std.fmt.allocPrintSentinel(allocator, "{s}.{s}", .{ fasta_path, case.suffix }, 0);
        try std.testing.expectEqual(std.os.linux.E.SUCCESS, std.os.linux.errno(std.os.linux.mknodat(std.os.linux.AT.FDCWD, sidecar, std.os.linux.S.IFIFO | 0o600, 0)));
        defer std.Io.Dir.cwd().deleteFile(io, sidecar) catch {};
        for ([_]main.index_format.LoadMode{ .lookup_full_map, .catalog, .stats_scan, .on_access }) |mode| {
            try expectLoadFailure(std.testing.allocator, fasta_path, mode, case.code);
        }
        try support.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "validate", "--index", fasta_path }, case.code);
    }
}

fn expectStreamMatchesCatalog(allocator: std.mem.Allocator, fasta_path: []const u8, selection: main.index_format.IndexSelection) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    var catalog = try main.index_format.loadIndex(allocator, io, fasta_path, selection, .catalog, &failure);
    defer catalog.deinit();
    const stream = try allocator.create(main.index_format.EntryStream);
    defer allocator.destroy(stream);
    try stream.open(allocator, io, fasta_path, selection, &failure);
    defer stream.close();
    for (0..catalog.records.len()) |i| {
        const record = try recordAt(&catalog, i);
        const entry = (try stream.next()).?;
        try std.testing.expectEqual(@as(u64, i + 1), entry.ordinal);
        try std.testing.expectEqualStrings(catalog.recordName(i).?, entry.name);
        try std.testing.expectEqualSlices(u8, std.mem.asBytes(&record), std.mem.asBytes(&entry.record));
        if (record.isUniformWidth()) continue;
        const lines = catalog.sideTableLines(record);
        try std.testing.expectEqual(@as(u64, lines.len), try stream.sideLineCount(entry.record));
        for (lines, 0..) |line, k| try std.testing.expectEqual(line, try stream.sideLine(entry.record, k));
    }
    try std.testing.expectEqual(@as(?*const main.index_format.EntryStream.Entry, null), try stream.next());
}

test "[property] - [index loader]: an entry stream yields the catalog's entries in order" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var prng = std.Random.DefaultPrng.init(0x1d5e_57ea);
    const random = prng.random();
    for (0..40) |_| {
        var data: std.ArrayList(u8) = .empty;
        const uniform = random.boolean();
        for (0..random.intRangeAtMost(usize, 1, 8)) |_| {
            try data.print(allocator, ">{s} d\n", .{([_][]const u8{ "a", "b", "chr1", "x" ** 65_000 })[random.uintLessThan(usize, if (random.uintLessThan(u8, 20) == 0) 4 else 3)]});
            const width = random.intRangeAtMost(usize, 1, 9);
            for (0..random.intRangeAtMost(usize, 1, 6)) |_| {
                const n = if (uniform) width else random.intRangeAtMost(usize, 1, 9);
                try data.appendNTimes(allocator, 'A', n);
                try data.appendSlice(allocator, if (random.uintLessThan(u8, 8) == 0 and !uniform) "\r\n" else "\n");
            }
        }
        const path = try writeFastaArtifact(allocator, "entry-stream", data.items);
        defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
        const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{path});
        defer std.Io.Dir.cwd().deleteFile(io, zfi_path) catch {};
        try support.writeZfi(allocator, path, data.items);
        try expectStreamMatchesCatalog(allocator, path, .{ .zfi = zfi_path });
        var fai: std.Io.Writer.Allocating = .init(allocator);
        _ = scanFaiData(allocator, data.items, &fai.writer) catch continue;
        const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{path});
        defer std.Io.Dir.cwd().deleteFile(io, fai_path) catch {};
        {
            const file = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
            defer file.close(io);
            try std.Io.File.writeStreamingAll(file, io, fai.written());
        }
        try expectStreamMatchesCatalog(allocator, path, .{ .fai = fai_path });
    }
}

test "[property] - [zfi name table]: every lookup returns the first record with that name" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var prng = std.Random.DefaultPrng.init(0x5a46_4954);
    const random = prng.random();
    for ([_]usize{ 1, 2, 3, 4, 7, 8, 9, 63, 64, 65, 1000 }) |count| {
        var fasta: std.ArrayList(u8) = .empty;
        var first: std.StringHashMapUnmanaged(usize) = .empty;
        for (0..count) |i| {
            var buffer: [32]u8 = undefined;
            const name = switch (random.uintLessThan(u8, 4)) {
                0 => "",
                1 => try std.fmt.bufPrint(&buffer, "dup{d}", .{random.uintLessThan(u8, 5)}),
                else => try std.fmt.bufPrint(&buffer, "name-{d}", .{i}),
            };
            const owned = try allocator.dupe(u8, name);
            const entry = try first.getOrPut(allocator, owned);
            if (!entry.found_existing) entry.value_ptr.* = i;
            try fasta.print(allocator, ">{s}\nACGT\n", .{owned});
        }
        var index = try main.indexer.scanZfiData(allocator, fasta.items);
        const bytes = try zfiBytes(allocator, &index, fasta.items.len, 0);
        const paths = try writeFastaAndRawZfi(allocator, "zfi-name-table", fasta.items, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        for ([_]main.index_format.LoadMode{ .on_access, .{ .requested = &.{} } }) |mode| {
            var loaded = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, mode);
            defer loaded.deinit();
            var failure: diagnostics.Failure = undefined;
            failure.reset(.get);
            var names = first.iterator();
            while (names.next()) |entry| {
                try std.testing.expectEqual(@as(?usize, entry.value_ptr.*), try loaded.findName(entry.key_ptr.*, &failure));
            }
            try std.testing.expectEqual(@as(?usize, null), try loaded.findName("absent", &failure));
        }
    }
}

test "[failure] - [zfi name table]: a damaged table fails in every mode that reads it" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n>b\nGG\n>c\nT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    const original = try zfiBytes(allocator, &index, fasta.len, 0);
    const table_len = 4 * main.index_format.NameTable.slotCount(3);
    const table_start = original.len - 24 - table_len;
    for ([_]enum { out_of_range, wrong_record }{ .out_of_range, .wrong_record }) |defect| {
        const bytes = try allocator.dupe(u8, original);
        var slot: usize = 0;
        while (std.mem.readInt(u32, bytes[table_start + slot * 4 ..][0..4], .little) == 0) slot += 1;
        const value = std.mem.readInt(u32, bytes[table_start + slot * 4 ..][0..4], .little);
        const ordinal_mask: u32 = 0b11;
        std.mem.writeInt(u32, bytes[table_start + slot * 4 ..][0..4], switch (defect) {
            .out_of_range => value & ~ordinal_mask,
            .wrong_record => (value & ~ordinal_mask) | (((value & ordinal_mask) % 3) + 1),
        }, .little);
        const paths = try writeFastaAndRawZfi(allocator, "zfi-damaged-table", fasta, bytes);
        defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
        defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
        for ([_]main.index_format.LoadMode{.catalog}) |mode| {
            var failure: diagnostics.Failure = undefined;
            try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.name_table_invalid", &failure);
            try std.testing.expectEqual(@as(u64, slot), failure.primary().?.details.@"index.name_table_invalid".slot);
            try std.testing.expectEqual(@as(u64, table_start + slot * 4), failure.primary().?.at.location.offset.?);
        }
        var stream: main.index_format.EntryStream = undefined;
        var failure: diagnostics.Failure = undefined;
        failure.reset(.validate);
        try stream.open(allocator, io, paths.fasta_path, .automatic, &failure);
        defer stream.close();
        const err = while (true) {
            _ = (stream.next() catch |err| break err) orelse return error.ExpectedIndexFailure;
        };
        try std.testing.expectEqual(error.InvalidIndex, err);
        try std.testing.expectEqual(diagnostics.Code.@"index.name_table_invalid", failure.primary().?.code());
        if (defect == .out_of_range) {
            var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .lookup_full_map);
            defer lazy.deinit();
            failure.reset(.get);
            var found_damage = false;
            for ([_][]const u8{ "a", "b", "c" }) |name| {
                _ = lazy.findName(name, &failure) catch {
                    try std.testing.expectEqual(diagnostics.Code.@"index.name_table_invalid", failure.primary().?.code());
                    found_damage = true;
                    continue;
                };
            }
            try std.testing.expect(found_damage);
        }
    }
}

test "[failure] - [zfi name table]: a file too short for the table is truncated" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    const original = try zfiBytes(allocator, &index, fasta.len, 0);
    var bytes: std.ArrayList(u8) = .empty;
    try bytes.appendSlice(allocator, original[0 .. 16 + 24]);
    try bytes.appendSlice(allocator, original[original.len - 24 ..]);
    std.mem.writeInt(u64, bytes.items[bytes.items.len - 8 ..][0..8], 0, .little);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-short-table", fasta, bytes.items);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    for ([_]main.index_format.LoadMode{ .on_access, .{ .requested = &.{"a"} } }) |mode| {
        var failure: diagnostics.Failure = undefined;
        try expectIndexLoadCode(allocator, paths.fasta_path, mode, .@"index.truncated", &failure);
        const facts = failure.primary().?.details.@"index.truncated";
        try std.testing.expectEqual(diagnostics.Section.name_table, facts.section);
        try std.testing.expectEqual(@as(u128, 40 + 4 * main.index_format.NameTable.slotCount(1)), facts.required_end);
        try std.testing.expectEqual(@as(u64, 40), facts.available_end);
    }
}

test "[failure] - [zfi name table]: on-demand lookups report damaged names and failed reads" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n>b\nGG\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    index.records.items[1].name_offset = 1000;
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-on-demand-faults", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var loaded = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{"b"} });
    defer loaded.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectEqual(@as(?usize, 0), try loaded.findName("a", &failure));
    try std.testing.expectError(error.InvalidIndex, loaded.findName("b", &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"index.name_out_of_bounds", primary.code());
    try std.testing.expectEqual(@as(?u64, 2), primary.at.record);
    try std.testing.expectEqual(@as(?u64, 16 + 24), primary.at.location.offset);

    const Fault = struct {
        var reads: usize = 0;
        fn read(userdata: ?*anyopaque, file: std.Io.File, data: []const []u8, offset: u64) std.Io.File.ReadPositionalError!usize {
            reads += 1;
            if (reads > 2) return error.InputOutput;
            return std.testing.io.vtable.fileReadPositional(userdata, file, data, offset);
        }
    };
    var vtable = io.vtable.*;
    vtable.fileReadPositional = Fault.read;
    const faulty_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    failure.reset(.get);
    var faulty = try main.index_format.loadIndex(allocator, faulty_io, paths.fasta_path, .automatic, .{ .requested = &.{"a"} }, &failure);
    defer faulty.deinit();
    try std.testing.expectError(error.Io, faulty.findName("a", &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.SourceId.zfi, failure.primary().?.at.source.?);
    try std.testing.expectEqual(@as(anyerror, error.InputOutput), failure.primary().?.cause.?.zig);
}

test "[property] - [zfi name table]: batch prefetch agrees with one lookup at a time" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var prng = std.Random.DefaultPrng.init(0x7072_6566);
    const random = prng.random();
    var fasta: std.ArrayList(u8) = .empty;
    var requests: std.ArrayList([]const u8) = .empty;
    for (0..3000) |i| {
        var buffer: [32]u8 = undefined;
        const name = try allocator.dupe(u8, switch (random.uintLessThan(u8, 5)) {
            0 => "",
            1 => try std.fmt.bufPrint(&buffer, "dup{d}", .{random.uintLessThan(u8, 40)}),
            else => try std.fmt.bufPrint(&buffer, "name-{d}", .{i}),
        });
        try fasta.print(allocator, ">{s}\n{s}", .{ name, if (i % 5 == 0) "ACG\nTACGT\nA\n" else "ACGT\nACGT\nA\n" });
        if (random.boolean()) try requests.append(allocator, name);
    }
    for (0..200) |i| try requests.append(allocator, try std.fmt.allocPrint(allocator, "absent-{d}", .{i}));
    random.shuffle([]const u8, requests.items);
    var index = try main.indexer.scanZfiData(allocator, fasta.items);
    const bytes = try zfiBytes(allocator, &index, fasta.items.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-prefetch", fasta.items, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};

    var full = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .on_access);
    defer full.deinit();
    var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{} });
    defer lazy.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    var buffer: [64 * 1024]u8 = undefined;
    lazy.prefetchNames(requests.items, &buffer);
    try std.testing.expectEqual(requests.items.len, lazy.on_demand.?.prefetched.len);
    for (lazy.on_demand.?.prefetched) |result| try std.testing.expect(result.state != .unknown);
    for (requests.items, 0..) |name, i| {
        const expected = full.lookupName(name);
        try std.testing.expectEqual(expected, try lazy.findNameAt(i, name, &failure));
        if (expected) |ordinal| {
            const rec = try lazy.checkedRecord(ordinal, &failure);
            try std.testing.expectEqual(try recordAt(&full, ordinal), rec);
            try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(full.sideTableLines(rec)), std.mem.sliceAsBytes(lazy.sideTableLines(rec)));
        }
    }
    lazy.clearPrefetch();
    try std.testing.expectEqual(@as(usize, 0), lazy.on_demand.?.prefetched.len);
}

test "[failure] - [zfi name table]: batch prefetch leaves a damaged side table to the record check" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACG\nTACGT\nA\n>b\nGG\nGGG\nG\n>c\nT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    const second = index.records.items[1].sideTableOffset() - (16 + index.records.items.len * 24);
    std.mem.writeInt(u32, index.side_tables.items[@intCast(second + 8 + 16)..][0..4], 0, .little);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-prefetch-sides", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{} });
    defer lazy.deinit();
    const names = [_][]const u8{ "a", "b", "c" };
    var buffer: [64 * 1024]u8 = undefined;
    lazy.prefetchNames(&names, &buffer);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectEqual(@as(?usize, 0), try lazy.findNameAt(0, names[0], &failure));
    try std.testing.expectEqual(index.records.items[0].seq_len, (try lazy.checkedRecord(0, &failure)).seq_len);
    try std.testing.expectEqual(@as(?usize, 1), try lazy.findNameAt(1, names[1], &failure));
    try std.testing.expectError(error.InvalidIndex, lazy.checkedRecord(1, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"index.side_line_geometry_invalid", failure.primary().?.code());
}

test "[edge] - [zfi name table]: batch prefetch hands a hash-bit match with another name to the sequential lookup" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n>b\nGG\n>c\nT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    // An entry at b's home slot with b's hash bits that points at a, as a colliding name's would.
    const slots = main.index_format.NameTable.slotCount(3);
    const table = bytes[bytes.len - 24 - slots * 4 ..][0 .. slots * 4];
    @memset(table, 0);
    const Entry = struct { name: []const u8, ordinal: u32 };
    for ([_]Entry{ .{ .name = "b", .ordinal = 0 }, .{ .name = "a", .ordinal = 0 }, .{ .name = "b", .ordinal = 1 }, .{ .name = "c", .ordinal = 2 } }) |entry| {
        const name_hash = main.index_format.NameTable.hash(entry.name);
        var slot: usize = @intCast((@as(u128, name_hash) * slots) >> 64);
        while (std.mem.readInt(u32, table[slot * 4 ..][0..4], .little) != 0) slot = (slot + 1) % slots;
        const ordinal_bits = 2;
        const value: u32 = @intCast((name_hash & ((1 << (32 - ordinal_bits)) - 1)) << ordinal_bits | (entry.ordinal + 1));
        std.mem.writeInt(u32, table[slot * 4 ..][0..4], value, .little);
    }
    const paths = try writeFastaAndRawZfi(allocator, "zfi-prefetch-collision", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{} });
    defer lazy.deinit();
    const names = [_][]const u8{ "a", "b", "c" };
    var buffer: [64 * 1024]u8 = undefined;
    lazy.prefetchNames(&names, &buffer);
    try std.testing.expect(lazy.on_demand.?.prefetched[1].state == .unknown);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    for (names, 0..) |name, i| try std.testing.expectEqual(@as(?usize, i), try lazy.findNameAt(i, name, &failure));
}

test "[failure] - [zfi name table]: batch prefetch leaves damage to the lookup that reaches it" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = ">a\nACGT\n>b\nGG\n>c\nT\n";
    var index = try main.indexer.scanZfiData(allocator, fasta);
    index.records.items[1].name_offset = 1000;
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-prefetch-damage", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{} });
    defer lazy.deinit();
    const names = [_][]const u8{ "a", "b", "c", "x" };
    var buffer: [64 * 1024]u8 = undefined;
    lazy.prefetchNames(&names, &buffer);
    try std.testing.expect(lazy.on_demand.?.prefetched[1].state == .unknown);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectEqual(@as(?usize, 0), try lazy.findNameAt(0, names[0], &failure));
    try std.testing.expectEqual(@as(?usize, 2), try lazy.findNameAt(2, names[2], &failure));
    try std.testing.expectEqual(@as(?usize, null), try lazy.findNameAt(3, names[3], &failure));
    try std.testing.expectError(error.InvalidIndex, lazy.findNameAt(1, names[1], &failure));
    try std.testing.expectEqual(diagnostics.Code.@"index.name_out_of_bounds", failure.primary().?.code());
}

test "[edge] - [zfi name table]: on-demand lookups compare names longer than one read" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const long = try allocator.alloc(u8, main.indexer.MAX_INDEX_NAME_LEN);
    @memset(long, 'n');
    const near = try allocator.dupe(u8, long);
    near[near.len - 1] = 'm';
    const fasta = try std.fmt.allocPrint(allocator, ">{s}\nACGT\n>{s}\nGG\n", .{ near, long });
    var index = try main.indexer.scanZfiData(allocator, fasta);
    const bytes = try zfiBytes(allocator, &index, fasta.len, 0);
    const paths = try writeFastaAndRawZfi(allocator, "zfi-long-lookup", fasta, bytes);
    defer std.Io.Dir.cwd().deleteFile(io, paths.fasta_path) catch {};
    defer std.Io.Dir.cwd().deleteFile(io, paths.zfi_path) catch {};
    var lazy = try support.loadIndex(std.testing.allocator, io, paths.fasta_path, .{ .requested = &.{} });
    defer lazy.deinit();
    var failure: diagnostics.Failure = undefined;
    failure.reset(.get);
    try std.testing.expectEqual(@as(?usize, 1), try lazy.findName(long, &failure));
    try std.testing.expectEqual(@as(?usize, 0), try lazy.findName(near, &failure));
    try std.testing.expectEqual(@as(?usize, null), try lazy.findName(long[1..], &failure));
}
