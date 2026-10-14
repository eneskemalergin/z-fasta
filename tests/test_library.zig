//! Public package tests: use only `@import("z-fasta")`, as an external consumer does.

const std = @import("std");

const zfasta = @import("z-fasta");

const diagnostics = zfasta.diagnostics;
const testing = std.testing;
const io = testing.io;

const ZFASTA_BIN = @import("utility.zig").ZFASTA_BIN;

const UNIFORM = ">a desc\nACGTA\nCGTAC\nGT\n>b\nTTTTT\nTT\n";
const NONUNIFORM = ">a desc\nACG\nTACGT\nACGT\n>b\nTTT\nTTTT\n";

const Fixture = struct {
    tmp: testing.TmpDir,
    root: [:0]u8,

    fn init(allocator: std.mem.Allocator) !Fixture {
        var tmp = testing.tmpDir(.{});
        errdefer tmp.cleanup();
        const root = try tmp.dir.realPathFileAlloc(io, ".", allocator);
        return .{ .tmp = tmp, .root = root };
    }

    fn deinit(self: *Fixture, allocator: std.mem.Allocator) void {
        allocator.free(self.root);
        self.tmp.cleanup();
    }

    fn path(self: *const Fixture, allocator: std.mem.Allocator, name: []const u8) ![]u8 {
        return std.fs.path.join(allocator, &.{ self.root, name });
    }

    fn write(self: *const Fixture, name: []const u8, bytes: []const u8) !void {
        try self.tmp.dir.writeFile(io, .{ .sub_path = name, .data = bytes });
    }
};

fn runCli(allocator: std.mem.Allocator, argv: []const []const u8, stdout_path: ?[]const u8) !void {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = argv, .stdout_limit = .limited(1 << 20) });
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    try testing.expectEqual(@as(u8, 0), result.term.exited);
    if (stdout_path) |out| try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = out, .data = result.stdout });
}

fn fetch(reader: *zfasta.IndexedReader, name: []const u8, interval: zfasta.Interval, dst: []u8) ![]u8 {
    var failure: diagnostics.Failure = undefined;
    const info = try reader.lookup(name, &failure) orelse return error.TestUnexpectedResult;
    return reader.fetchInto(info.ordinal, interval, dst, &failure);
}

fn expectCode(code: diagnostics.Code, failure: *const diagnostics.Failure) !void {
    try testing.expectEqual(code, failure.primary().?.code());
}

fn expectRecordsMatchIndex(allocator: std.mem.Allocator, fixture: *const Fixture, name: []const u8, bytes: []const u8) !void {
    try fixture.write(name, bytes);
    const fasta = try fixture.path(allocator, name);
    defer allocator.free(fasta);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", fasta }, null);

    var index_failure: diagnostics.Failure = undefined;
    var indexed = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &index_failure);
    defer indexed.deinit();
    var source: std.Io.Reader = .fixed(bytes);
    var storage: [256]u8 = undefined;
    var failure: diagnostics.Failure = undefined;
    var records = zfasta.RecordReader.init(&source, &storage, &failure);
    var ordinal: usize = 0;
    while (try records.next()) |record| : (ordinal += 1) {
        const info = try indexed.record(ordinal, &index_failure) orelse return error.TestUnexpectedResult;
        try testing.expectEqualStrings(info.name, record.id);
        var bases: [256]u8 = undefined;
        try testing.expectEqualStrings(try indexed.fetchInto(ordinal, .{ .start = 0, .end = info.length }, &bases, &index_failure), record.sequence);
    }
    try testing.expectEqual(indexed.recordCount(), ordinal);
}

fn expectCopiedBases(fasta: []const u8, want: []const u8, storage: []u8, chunk: ?usize) !void {
    var whole: std.Io.Reader = .fixed(fasta);
    var backing: [64]u8 = undefined;
    const calls = [_]std.testing.Reader.Call{.{ .buffer = fasta }};
    var split = std.testing.Reader.init(&backing, &calls);
    if (chunk) |limit| split.artificial_limit = .limited(limit);
    var failure: diagnostics.Failure = undefined;
    var records = zfasta.RecordReader.init(if (chunk == null) &whole else &split.interface, storage, &failure);
    if (want.len == 0) return testing.expect((try records.next()) == null);
    if (storage.len < 1 + want.len) {
        try testing.expectError(error.CapacityExceeded, records.next());
        const observed = failure.primary().?.details.@"resource.line_buffer_exceeded".observed_bytes.at_least;
        try testing.expect(observed > storage.len and observed <= 1 + want.len);
        return;
    }
    const record = (try records.next()) orelse return error.TestUnexpectedResult;
    try testing.expectEqualStrings("x", record.id);
    try testing.expectEqualStrings(want, record.sequence);
    try testing.expect((try records.next()) == null);
}

fn reportCase(err: anyerror, case: usize, mode: usize, body_len: usize, split: ?usize) anyerror {
    std.debug.print("record reader mismatch: seed={d}, case={d}, mode={d}, body={d}, split={?d}\n", .{ std.testing.random_seed, case, mode, body_len, split });
    return err;
}

const FailingSource = struct {
    interface: std.Io.Reader,
    data: []const u8,

    fn init(data: []const u8, buffer: []u8) FailingSource {
        return .{ .interface = .{ .vtable = &.{ .stream = stream }, .buffer = buffer, .seek = 0, .end = 0 }, .data = data };
    }

    fn stream(reader: *std.Io.Reader, writer: *std.Io.Writer, limit: std.Io.Limit) std.Io.Reader.StreamError!usize {
        const self: *FailingSource = @alignCast(@fieldParentPtr("interface", reader));
        if (self.data.len == 0) return error.ReadFailed;
        const count = try writer.write(limit.sliceConst(self.data));
        self.data = self.data[count..];
        return count;
    }
};

test "[integration] - [indexed reader]: one open reader serves repeated fetches" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("ref.fa", UNIFORM);
    const fasta = try fixture.path(allocator, "ref.fa");
    defer allocator.free(fasta);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", fasta }, null);

    var failure: diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
    defer reader.deinit();
    try testing.expectEqual(@as(usize, 2), reader.recordCount());
    const a = (try reader.lookup("a", &failure)).?;
    try testing.expectEqualStrings("a", a.name);
    try testing.expectEqual(@as(u64, 12), a.length);

    var buf: [16]u8 = undefined;
    try testing.expectEqualStrings("ACGTA", try reader.fetchInto(a.ordinal, .{ .start = 4, .end = 9 }, &buf, &failure));
    try testing.expectEqualStrings("ACGTACGTACGT", try reader.fetchInto(a.ordinal, .{ .start = 0, .end = 12 }, &buf, &failure));
    try testing.expectEqualStrings("TTTTTTT", try fetch(&reader, "b", .{ .start = 0, .end = 7 }, &buf));
    try testing.expect(try reader.lookup("missing", &failure) == null);
    try testing.expect(try reader.record(2, &failure) == null);
}

test "[property] - [indexed reader]: FAI, uniform ZFI and non-uniform ZFI agree" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("uniform.fa", UNIFORM);
    try fixture.write("nonuniform.fa", NONUNIFORM);
    const uniform = try fixture.path(allocator, "uniform.fa");
    defer allocator.free(uniform);
    const uniform_fai = try fixture.path(allocator, "uniform.fa.fai");
    defer allocator.free(uniform_fai);
    const nonuniform = try fixture.path(allocator, "nonuniform.fa");
    defer allocator.free(nonuniform);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", uniform }, null);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", "--emit-fai", uniform }, uniform_fai);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", nonuniform }, null);

    const Case = struct { path: []const u8, selection: zfasta.IndexSelection };
    const cases = [_]Case{
        .{ .path = uniform, .selection = .{ .fai = uniform_fai } },
        .{ .path = uniform, .selection = .automatic },
        .{ .path = nonuniform, .selection = .automatic },
    };
    const expected = "ACGTACGTACGT";
    for (cases) |case| {
        var failure: diagnostics.Failure = undefined;
        var reader = try zfasta.IndexedReader.open(allocator, io, case.path, .{ .index = case.selection }, &failure);
        defer reader.deinit();
        var start: u64 = 0;
        while (start <= expected.len) : (start += 1) {
            var end = start;
            while (end <= expected.len) : (end += 1) {
                var buf: [16]u8 = undefined;
                try testing.expectEqualStrings(expected[start..end], try fetch(&reader, "a", .{ .start = start, .end = end }, &buf));
            }
        }
    }
}

test "[edge] - [indexed reader]: strict intervals fail before changing the destination" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("ref.fa", UNIFORM);
    const fasta = try fixture.path(allocator, "ref.fa");
    defer allocator.free(fasta);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", fasta }, null);

    var failure: diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
    defer reader.deinit();
    var buf: [4]u8 = "xxxx".*;
    try testing.expectEqualStrings("", try reader.fetchInto(0, .{ .start = 12, .end = 12 }, &buf, &failure));

    try testing.expectError(error.InvalidRequest, reader.fetchInto(2, .{ .start = 0, .end = 1 }, &buf, &failure));
    try expectCode(.@"request.ordinal_out_of_range", &failure);
    try testing.expectError(error.InvalidRequest, reader.fetchInto(0, .{ .start = 3, .end = 2 }, &buf, &failure));
    try expectCode(.@"request.end_before_start", &failure);
    try testing.expectError(error.InvalidRequest, reader.fetchInto(0, .{ .start = 0, .end = 13 }, &buf, &failure));
    try expectCode(.@"request.end_after_sequence", &failure);
    try testing.expectError(error.InvalidRequest, reader.fetchInto(0, .{ .start = 0, .end = 5 }, &buf, &failure));
    try expectCode(.@"request.output_buffer_too_small", &failure);
    try testing.expectEqual(@as(u64, 5), failure.primary().?.details.@"request.output_buffer_too_small".required_bytes);
    try testing.expectEqualStrings("xxxx", &buf);
}

test "[integration] - [indexed reader]: lookup is first-wins while ordinals reach duplicates" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("dup.fa", ">x\nAAAA\n>x\nCCCC\n");
    const fasta = try fixture.path(allocator, "dup.fa");
    defer allocator.free(fasta);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", fasta }, null);

    var failure: diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
    defer reader.deinit();
    try testing.expectEqual(@as(usize, 2), reader.recordCount());
    try testing.expectEqual(@as(usize, 0), (try reader.lookup("x", &failure)).?.ordinal);
    var buf: [4]u8 = undefined;
    try testing.expectEqualStrings("CCCC", try reader.fetchInto(1, .{ .start = 0, .end = 4 }, &buf, &failure));
}

test "[failure] - [indexed reader]: opens a damaged .zfi and reports the record it reads" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("ref.fa", UNIFORM);
    const fasta = try fixture.path(allocator, "ref.fa");
    defer allocator.free(fasta);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", fasta }, null);
    const zfi = try fixture.tmp.dir.readFileAlloc(io, "ref.fa.zfi", allocator, .unlimited);
    defer allocator.free(zfi);
    // The second record's name length, bits 48-63 of its first word, now reaches past the names.
    zfi[16 + 24 + 6] = 0xff;
    zfi[16 + 24 + 7] = 0xff;
    try fixture.write("ref.fa.zfi", zfi);

    var failure: diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
    defer reader.deinit();
    try testing.expectEqualStrings("a", (try reader.record(0, &failure)).?.name);
    try testing.expectError(error.InvalidIndex, reader.record(1, &failure));
    try expectCode(.@"index.name_out_of_bounds", &failure);
    var buf: [8]u8 = undefined;
    try testing.expectError(error.InvalidIndex, reader.fetchInto(1, .{ .start = 0, .end = 7 }, &buf, &failure));
    try testing.expectEqualStrings("ACGTA", try fetch(&reader, "a", .{ .start = 0, .end = 5 }, &buf));
}

test "[failure] - [indexed reader]: explicit selection never falls back" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("ref.fa", UNIFORM);
    try fixture.write("ref.fa.zfi", "not a zfi index");
    const fasta = try fixture.path(allocator, "ref.fa");
    defer allocator.free(fasta);
    const fai = try fixture.path(allocator, "ref.fa.fai");
    defer allocator.free(fai);
    const missing = try fixture.path(allocator, "absent.fai");
    defer allocator.free(missing);
    try runCli(allocator, &.{ ZFASTA_BIN, "index", "--emit-fai", fasta }, fai);

    var failure: diagnostics.Failure = undefined;
    try testing.expectError(error.InvalidIndex, zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure));
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{ .index = .{ .fai = fai } }, &failure);
    reader.deinit();
    try testing.expectError(error.Io, zfasta.IndexedReader.open(allocator, io, fasta, .{ .index = .{ .fai = missing } }, &failure));
    try expectCode(.@"io.file_not_found", &failure);
    try testing.expectEqual(diagnostics.SourceId.fai, failure.primary().?.at.source.?);
}

test "[failure] - [indexed reader]: rejects gzip sources and whitespace in dense rows" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    try fixture.write("gz.fa", "\x1f\x8b\x08\x00\x00\x00\x00\x00");
    try fixture.write("gz.fa.fai", "a\t4\t3\t4\t5\n");
    try fixture.write("ws.fa", ">a\nAC T\n");
    try fixture.write("ws.fa.fai", "a\t4\t3\t4\t5\n");
    const gz = try fixture.path(allocator, "gz.fa");
    defer allocator.free(gz);
    const ws = try fixture.path(allocator, "ws.fa");
    defer allocator.free(ws);

    var failure: diagnostics.Failure = undefined;
    try testing.expectError(error.InvalidFasta, zfasta.IndexedReader.open(allocator, io, gz, .{}, &failure));
    try expectCode(.@"input.unsupported_format", &failure);

    var reader = try zfasta.IndexedReader.open(allocator, io, ws, .{}, &failure);
    defer reader.deinit();
    var buf: [4]u8 = undefined;
    try testing.expectError(error.InvalidIndex, reader.fetchInto(0, .{ .start = 0, .end = 4 }, &buf, &failure));
    try expectCode(.@"pair.region_length_mismatch", &failure);
    try testing.expectEqual(@as(u64, 3), failure.primary().?.details.@"pair.region_length_mismatch".observed_symbols);
}

// --- Record reader ---

test "[integration] - [record reader]: agrees with index on records, IDs, and sequences" {
    const allocator = testing.allocator;
    var fixture = try Fixture.init(allocator);
    defer fixture.deinit(allocator);
    const cases = [_]struct { name: []const u8, bytes: []const u8 }{
        .{ .name = "uniform.fa", .bytes = UNIFORM },
        .{ .name = "nonuniform.fa", .bytes = NONUNIFORM },
        .{ .name = "crlf.fa", .bytes = ">a x\r\nACGT\r\nAC\r\n>b\r\nGG\r\n" },
        .{ .name = "empty.fa", .bytes = ">e1\n>a\nACGT\n>e2\n\n>b\nA\n>e3\n" },
        .{ .name = "dup.fa", .bytes = ">x\nAAAA\n>x\nCCCC\n" },
        .{ .name = "inner.fa", .bytes = ">a\nAC>GT\nA>\n>b\nG\n" },
        .{ .name = "noid.fa", .bytes = ">\nACGT\n> sp\nGG\n" },
        .{ .name = "spaces.fa", .bytes = ">a\nAC GT\nA\tC \n>b\nG" },
    };
    for (cases) |case| try expectRecordsMatchIndex(allocator, &fixture, case.name, case.bytes);
}

test "[property] - [record reader]: unbuffered and one-byte reads give the same records" {
    const input = ">a first desc\r\nAC GT\r\n  \r\nAC\tGT\n>b\n>c x\ty \r\nG>T\n\n>d\r\nAAAA";
    const expected = [_][3][]const u8{
        .{ "a", "first desc", "ACGTACGT" },
        .{ "c", "x\ty", "G>T" },
        .{ "d", "", "AAAA" },
    };
    for ([_]usize{ 0, 1, 2, 64 }) |capacity| {
        for ([_]usize{ 1, 3, 7 }) |chunk| {
            var backing: [64]u8 = undefined;
            const calls = [_]std.testing.Reader.Call{.{ .buffer = input }};
            var split = std.testing.Reader.init(backing[0..capacity], &calls);
            split.artificial_limit = .limited(chunk);
            var storage: [64]u8 = undefined;
            var failure: diagnostics.Failure = undefined;
            var records = zfasta.RecordReader.init(&split.interface, &storage, &failure);
            for (expected) |want| {
                const got = (try records.next()) orelse return error.TestUnexpectedResult;
                try testing.expectEqualStrings(want[0], got.id);
                try testing.expectEqualStrings(want[1], got.description);
                try testing.expectEqualStrings(want[2], got.sequence);
            }
            try testing.expect((try records.next()) == null);
        }
    }
}

test "[property] - [record reader]: copied bases equal a scalar filter for every size and split" {
    var prng = std.Random.DefaultPrng.init(std.testing.random_seed);
    const random = prng.random();
    var input: [3 + 1024]u8 = undefined;
    var expected: [1024]u8 = undefined;
    var storage: [1 + 1024]u8 = undefined;
    @memcpy(input[0..3], ">x\n");
    for (0..400) |case| {
        const mode = case % 4;
        const body_len = random.uintAtMost(usize, 1024);
        var count: usize = 0;
        for (input[3..][0..body_len], 0..) |*byte, i| {
            byte.* = switch (mode) {
                0 => "ACGTN"[random.uintLessThan(usize, 5)],
                1 => if (i % 61 == 60) '\n' else "ACGT"[random.uintLessThan(usize, 4)],
                2 => "ACGT \t\r\n"[random.uintLessThan(usize, 8)],
                else => arbitrary: {
                    const value = random.int(u8);
                    break :arbitrary if (value == '>') 'A' else value;
                },
            };
            if (byte.* > ' ') {
                expected[count] = byte.*;
                count += 1;
            }
        }
        const fasta = input[0 .. 3 + body_len];
        const want = expected[0..count];
        const split = 1 + random.uintLessThan(usize, 97);
        expectCopiedBases(fasta, want, storage[0 .. 1 + count], null) catch |err| return reportCase(err, case, mode, body_len, null);
        expectCopiedBases(fasta, want, storage[0 .. 1 + count], split) catch |err| return reportCase(err, case, mode, body_len, split);
        if (count > 0) expectCopiedBases(fasta, want, storage[0..count], null) catch |err| return reportCase(err, case, mode, body_len, null);
    }
}

test "[unit] - [record reader]: splits headers like index and trims descriptions" {
    const input = ">a  desc  text \r\nA\n>b\tx\nA\n> lead\nA\n>c\rmid\nA\n>\nA\n";
    const expected = [_][2][]const u8{ .{ "a", "desc  text" }, .{ "b", "x" }, .{ "", "lead" }, .{ "c", "mid" }, .{ "", "" } };
    var source: std.Io.Reader = .fixed(input);
    var storage: [32]u8 = undefined;
    var failure: diagnostics.Failure = undefined;
    var records = zfasta.RecordReader.init(&source, &storage, &failure);
    for (expected) |want| {
        const got = (try records.next()) orelse return error.TestUnexpectedResult;
        try testing.expectEqualStrings(want[0], got.id);
        try testing.expectEqualStrings(want[1], got.description);
        try testing.expectEqualStrings("A", got.sequence);
    }
    try testing.expect((try records.next()) == null);
}

test "[edge] - [record reader]: inputs without bases yield no records and complete" {
    for ([_][]const u8{ "", ">a\n", ">a\n>b\n\n \t\r\n" }) |input| {
        var source: std.Io.Reader = .fixed(input);
        var storage: [16]u8 = undefined;
        var failure: diagnostics.Failure = undefined;
        var records = zfasta.RecordReader.init(&source, &storage, &failure);
        try testing.expect((try records.next()) == null);
        try testing.expect((try records.next()) == null);
        try testing.expect(failure.complete);
    }
}

test "[unit] - [record reader]: reads UniProt-shaped protein records" {
    const line1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK";
    const line2 = "ALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWE";
    const line3 = "RVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDL";
    const input = ">sp|P99999|TEST_HUMAN Test protein one OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=2\n" ++
        line1 ++ "\n" ++ line2 ++ "\n" ++ line3 ++ "\n" ++
        ">tr|A0A000|A0A000_9ZZZZ Uncharacterized protein OS=Example organism OX=1\nMNNNNNNNNNKLLX*\n";
    var source: std.Io.Reader = .fixed(input);
    var storage: [512]u8 = undefined;
    var failure: diagnostics.Failure = undefined;
    var records = zfasta.RecordReader.init(&source, &storage, &failure);
    const first = (try records.next()).?;
    try testing.expectEqualStrings("sp|P99999|TEST_HUMAN", first.id);
    try testing.expectEqualStrings("Test protein one OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=2", first.description);
    try testing.expectEqualStrings(line1 ++ line2 ++ line3, first.sequence);
    const second = (try records.next()).?;
    try testing.expectEqualStrings("tr|A0A000|A0A000_9ZZZZ", second.id);
    try testing.expectEqualStrings("Uncharacterized protein OS=Example organism OX=1", second.description);
    try testing.expectEqualStrings("MNNNNNNNNNKLLX*", second.sequence);
    try testing.expect((try records.next()) == null);
}

test "[failure] - [record reader]: content before the first header fails at offset zero" {
    for ([_][]const u8{ "ACGT\n>a\nA\n", "\xEF\xBB\xBF>a\nA\n", "\n>a\nA\n" }) |input| {
        var source: std.Io.Reader = .fixed(input);
        var storage: [16]u8 = undefined;
        var failure: diagnostics.Failure = undefined;
        var records = zfasta.RecordReader.init(&source, &storage, &failure);
        try testing.expectError(error.InvalidFasta, records.next());
        try expectCode(.@"fasta.content_before_header", &failure);
        try testing.expectEqual(@as(?u64, 0), failure.primary().?.at.location.offset);
        try testing.expectError(error.InvalidFasta, records.next());
        try testing.expectEqual(@as(u8, 0), failure.secondary_len);
    }
}

test "[failure] - [record reader]: a record larger than storage fails and stays failed" {
    const Case = struct { input: []const u8, observed: u128 };
    for ([_]Case{ .{ .input = ">abcdef\nA\n", .observed = 6 }, .{ .input = ">a\nACGTACGT\n", .observed = 9 } }) |case| {
        var source: std.Io.Reader = .fixed(case.input);
        var storage: [4]u8 = undefined;
        var failure: diagnostics.Failure = undefined;
        var records = zfasta.RecordReader.init(&source, &storage, &failure);
        try testing.expectError(error.CapacityExceeded, records.next());
        try expectCode(.@"resource.line_buffer_exceeded", &failure);
        const facts = failure.primary().?.details.@"resource.line_buffer_exceeded";
        try testing.expectEqual(@as(u64, 4), facts.buffer_bytes);
        try testing.expectEqual(diagnostics.Measurement{ .at_least = case.observed }, facts.observed_bytes);
        try testing.expectEqual(@as(?u64, 1), failure.primary().?.at.record);
        try testing.expectError(error.CapacityExceeded, records.next());
        try testing.expectEqual(@as(u8, 0), failure.secondary_len);
    }
}

test "[failure] - [record reader]: a failing source reports the record and offset" {
    for ([_]usize{ 0, 8 }) |capacity| {
        var buffer: [8]u8 = undefined;
        var failing = FailingSource.init(">a\nAC", buffer[0..capacity]);
        var storage: [16]u8 = undefined;
        var failure: diagnostics.Failure = undefined;
        var records = zfasta.RecordReader.init(&failing.interface, &storage, &failure);
        try testing.expectError(error.Io, records.next());
        try expectCode(.@"io.read_failed", &failure);
        try testing.expectEqual(@as(?u64, 1), failure.primary().?.at.record);
        try testing.expectEqual(@as(?u64, 5), failure.primary().?.at.location.offset);
    }
}

// --- Validation ---

test "[unit] - [validation]: validates and repairs a stream through the package" {
    const data = ">a desc\r\nAC GT\r\nACG\n>a\nAC?\n";
    const report = try testing.allocator.create(zfasta.validation.Report);
    defer testing.allocator.destroy(report);
    var buffer: [16]u8 = undefined;
    var failure: diagnostics.Failure = undefined;

    var reader = std.Io.Reader.fixed(data);
    try zfasta.validation.validate(testing.allocator, &reader, &buffer, null, .{}, report, &failure);
    try testing.expect(failure.complete);
    try testing.expectEqual(@as(u64, 2), report.examined.records);
    try testing.expectEqual(@as(u64, 1), report.finding(.@"fasta.duplicate_identifier").occurrences);
    try testing.expectEqual(@as(u64, 1), report.finding(.@"sequence.invalid_symbol").occurrences);
    try testing.expectEqual(@as(u8, 1), report.status(false, false));

    reader = std.Io.Reader.fixed(data);
    var output: std.Io.Writer.Allocating = .init(testing.allocator);
    defer output.deinit();
    try zfasta.validation.repair(testing.allocator, &reader, &buffer, null, .{}, 4, &output.writer, report, &failure);
    try testing.expectEqualStrings(">a desc\nACGT\nACG\n>a\nAC?\n", output.written());

    reader = std.Io.Reader.fixed(">a\x01\nAC\n");
    output.clearRetainingCapacity();
    try testing.expectError(error.RepairRefused, zfasta.validation.repair(testing.allocator, &reader, &buffer, null, .{}, 60, &output.writer, report, &failure));
    try expectCode(.@"repair.blocked", &failure);
}
