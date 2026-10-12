//! Copies declared fixtures and prepares their indexes in the build's test directory.

const std = @import("std");
const commands = @import("main");

pub fn main(init: std.process.Init.Minimal) !void {
    const allocator = std.heap.page_allocator;
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    var args = try std.process.Args.Iterator.initAllocator(init.args, allocator);
    defer args.deinit();
    _ = args.skip();
    const root = args.next() orelse return error.MissingOutputDirectory;
    var directory = try std.Io.Dir.cwd().createDirPathOpen(io, root, .{});
    defer directory.close(io);
    try directory.createDirPath(io, "tests/data/gates");
    while (args.next()) |destination| {
        const source = args.next() orelse return error.MissingFixtureSource;
        try std.Io.Dir.cwd().copyFile(source, directory, destination, io, .{});
    }
    for ([_][]const u8{ "simple.fasta", "proteome.fasta", "edge_cases.fasta" }) |name| {
        const path = try std.fmt.allocPrint(allocator, "{s}/tests/data/{s}", .{ root, name });
        defer allocator.free(path);
        var failure: commands.index_format.diagnostics.Failure = undefined;
        _ = commands.indexer.runIndex(std.heap.page_allocator, io, init.environ, path, .{}, &failure) catch |err| {
            var buffer: [4096]u8 = undefined;
            var output = std.Io.File.stderr().writer(io, &buffer);
            try commands.reporting.writeText(&output.interface, &failure);
            try output.flush();
            return err;
        };
    }
}
