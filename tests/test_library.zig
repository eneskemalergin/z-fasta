//! Public package tests: use only `@import("z-fasta")`, as an external consumer does.

const std = @import("std");

const test_options = @import("test_options");
const zfasta = @import("z-fasta");

const diagnostics = zfasta.diagnostics;
const testing = std.testing;
const io = testing.io;

const ZFASTA_BIN = if (std.fs.path.isAbsolute(test_options.binary_path)) test_options.binary_path else test_options.build_root ++ std.fs.path.sep_str ++ test_options.binary_path;

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
    const info = reader.lookup(name) orelse return error.TestUnexpectedResult;
    return reader.fetchInto(info.ordinal, interval, dst, &failure);
}

fn expectCode(code: diagnostics.Code, failure: *const diagnostics.Failure) !void {
    try testing.expectEqual(code, failure.primary().?.code());
}

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
    const a = reader.lookup("a").?;
    try testing.expectEqualStrings("a", a.name);
    try testing.expectEqual(@as(u64, 12), a.length);

    var buf: [16]u8 = undefined;
    try testing.expectEqualStrings("ACGTA", try reader.fetchInto(a.ordinal, .{ .start = 4, .end = 9 }, &buf, &failure));
    try testing.expectEqualStrings("ACGTACGTACGT", try reader.fetchInto(a.ordinal, .{ .start = 0, .end = 12 }, &buf, &failure));
    try testing.expectEqualStrings("TTTTTTT", try fetch(&reader, "b", .{ .start = 0, .end = 7 }, &buf));
    try testing.expect(reader.lookup("missing") == null);
    try testing.expect(reader.record(2) == null);
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
    try runCli(allocator, &.{ ZFASTA_BIN, "index", "--no-dedup", fasta }, null);

    var failure: diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta, .{}, &failure);
    defer reader.deinit();
    try testing.expectEqual(@as(usize, 2), reader.recordCount());
    try testing.expectEqual(@as(usize, 0), reader.lookup("x").?.ordinal);
    var buf: [4]u8 = undefined;
    try testing.expectEqualStrings("CCCC", try reader.fetchInto(1, .{ .start = 0, .end = 4 }, &buf, &failure));
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
