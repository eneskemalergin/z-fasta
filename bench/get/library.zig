//! The package lane of the GET benchmark. `z-fasta-library FASTA REGIONS` opens the reference once
//! with `IndexedReader` and fetches every `NAME:START-END` line (one-based, inclusive) of REGIONS;
//! `z-fasta-library --records FASTA BYTES` streams every record with `RecordReader` through BYTES
//! of storage. Both print as `z-fasta get` prints, 60 symbols per line.

const std = @import("std");
const zfasta = @import("z-fasta");

pub fn main(init: std.process.Init.Minimal) !void {
    const allocator = std.heap.page_allocator;
    var args = try std.process.Args.Iterator.initAllocator(init.args, allocator);
    defer args.deinit();
    _ = args.skip();
    const first = args.next() orelse return error.MissingFasta;
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    var buffer: [64 * 1024]u8 = undefined;
    var stdout = std.Io.File.Writer.initStreaming(.stdout(), io, &buffer);
    if (std.mem.eql(u8, first, "--records")) {
        const fasta_path = args.next() orelse return error.MissingFasta;
        const storage_bytes = try std.fmt.parseInt(usize, args.next() orelse return error.MissingStorage, 10);
        try streamRecords(allocator, io, fasta_path, storage_bytes, &stdout.interface);
    } else {
        try fetchRegions(allocator, io, first, args.next() orelse return error.MissingRegions, &stdout.interface);
    }
    try stdout.interface.flush();
}

fn fetchRegions(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, regions_path: []const u8, out: *std.Io.Writer) !void {
    var failure: zfasta.diagnostics.Failure = undefined;
    var reader = try zfasta.IndexedReader.open(allocator, io, fasta_path, .{}, &failure);
    defer reader.deinit();
    const regions = try std.Io.Dir.cwd().readFileAlloc(io, regions_path, allocator, .unlimited);
    defer allocator.free(regions);
    var bases: std.ArrayList(u8) = .empty;
    defer bases.deinit(allocator);
    var lines = std.mem.tokenizeScalar(u8, regions, '\n');
    while (lines.next()) |line| {
        const colon = std.mem.lastIndexOfScalar(u8, line, ':') orelse return error.InvalidRegion;
        const dash = std.mem.indexOfScalarPos(u8, line, colon, '-') orelse return error.InvalidRegion;
        const start = try std.fmt.parseInt(u64, line[colon + 1 .. dash], 10);
        const end = try std.fmt.parseInt(u64, line[dash + 1 ..], 10);
        if (start == 0 or end < start) return error.InvalidRegion;
        const record = reader.lookup(line[0..colon]) orelse return error.UnknownRecord;
        try bases.resize(allocator, @intCast(end - start + 1));
        try writeRecord(out, line, try reader.fetchInto(record.ordinal, .{ .start = start - 1, .end = end }, bases.items, &failure));
    }
}

fn streamRecords(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, storage_bytes: usize, out: *std.Io.Writer) !void {
    const file = try std.Io.Dir.cwd().openFile(io, fasta_path, .{});
    defer file.close(io);
    var read_buffer: [256 * 1024]u8 = undefined;
    var source = file.reader(io, &read_buffer);
    const storage = try allocator.alloc(u8, storage_bytes);
    defer allocator.free(storage);
    var failure: zfasta.diagnostics.Failure = undefined;
    var records = zfasta.RecordReader.init(&source.interface, storage, &failure);
    while (try records.next()) |record| try writeRecord(out, record.id, record.sequence);
}

fn writeRecord(out: *std.Io.Writer, header: []const u8, sequence: []const u8) !void {
    try out.print(">{s}\n", .{header});
    var rest = sequence;
    while (rest.len != 0) {
        const take = @min(rest.len, 60);
        try out.print("{s}\n", .{rest[0..take]});
        rest = rest[take..];
    }
}
