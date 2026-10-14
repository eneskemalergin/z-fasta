//! Separate program that calls the four commands' entry points without the CLI dispatcher, first
//! with failing input and then with valid input.

const std = @import("std");
const commands = @import("main");

pub fn main(init: std.process.Init.Minimal) !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer if (debug_allocator.deinit() == .leak) @panic("leaked memory");
    const allocator = debug_allocator.allocator();
    var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
    defer threaded.deinit();
    backend = threaded.io();
    var vtable = backend.vtable.*;
    vtable.operate = discardStdout;
    const io = std.Io{ .userdata = backend.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    const missing = "tests/data/diagnostics-caller-does-not-exist.fa";
    const overlong: [std.Io.Dir.max_path_bytes + 1]u8 = @splat('x');
    const valid = "tests/data/simple.fasta";
    const get_options = commands.getter.Options{ .source = .{ .positional = &.{"seq1:1-4"} } };
    for ([_][]const u8{ missing, &overlong }) |path| {
        for ([_]bool{ false, true }) |emit_fai| {
            try std.testing.expectError(error.Io, commands.indexer.runIndex(allocator, io, init.environ, path, .{ .emit_fai = emit_fai }, &failure));
            try expectOpenFailure(&failure, .index, path);
        }
        try std.testing.expectError(error.Io, commands.getter.runGet(allocator, io, path, get_options, &failure));
        try expectOpenFailure(&failure, .get, path);
        try std.testing.expectError(error.Io, commands.stats.runStats(allocator, io, path, &failure));
        try expectOpenFailure(&failure, .stats, path);
        try std.testing.expectError(error.Io, commands.validator.runValidate(allocator, io, path, .{}, &failure));
        try expectOpenFailure(&failure, .validate, path);
    }

    const indexed = try commands.indexer.runIndex(allocator, io, init.environ, valid, .{ .emit_fai = true }, &failure);
    try std.testing.expectEqual(@as(u64, 2), indexed.record_count);
    try expectComplete(&failure);

    const extracted = try commands.getter.runGet(allocator, io, valid, get_options, &failure);
    try std.testing.expectEqual(@as(u64, 4), extracted.total_bases);
    try expectComplete(&failure);

    try commands.stats.runStats(allocator, io, valid, &failure);
    try expectComplete(&failure);

    const validated = try commands.validator.runValidate(allocator, io, valid, .{}, &failure);
    try std.testing.expectEqual(@as(u8, 0), validated.status);
    try expectComplete(&failure);
}

const diagnostics = commands.index_format.diagnostics;

var backend: std.Io = undefined;

fn discardStdout(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
    if (operation == .file_write_streaming) {
        const write = operation.file_write_streaming;
        if (write.file.handle == std.Io.File.stdout().handle) {
            var size = write.header.len;
            for (write.data, 0..) |data, i| size += data.len * (if (i + 1 == write.data.len) write.splat else 1);
            return .{ .file_write_streaming = size };
        }
    }
    return backend.vtable.operate(userdata, operation);
}

// Fails at open, before any output: `path` names no file, or is longer than the OS accepts.
fn expectOpenFailure(failure: *const diagnostics.Failure, operation: diagnostics.Operation, path: []const u8) !void {
    const primary = failure.primary().?;
    try std.testing.expectEqual(operation, failure.operation);
    try std.testing.expectEqual(diagnostics.Stage.open, primary.stage);
    try std.testing.expectEqual(diagnostics.SourceId.fasta, primary.at.source);
    try std.testing.expectEqualStrings(path, failure.sourceLabel(.fasta).?);
    try std.testing.expectEqual(diagnostics.SourceKind.path, failure.sources.get(.fasta).?.kind);
    if (path.len > std.Io.Dir.max_path_bytes) {
        try std.testing.expectEqual(diagnostics.Code.@"io.invalid_path", primary.code());
        try std.testing.expectEqual(error.NameTooLong, primary.cause.?.zig);
        const details = primary.details.@"io.invalid_path";
        try std.testing.expectEqual(diagnostics.PathReason.too_long, details.reason);
        try std.testing.expect(details.limit_bytes == null);
    } else {
        try std.testing.expectEqual(diagnostics.Code.@"io.file_not_found", primary.code());
        try std.testing.expectEqual(error.FileNotFound, primary.cause.?.zig);
    }
    try std.testing.expect(!failure.complete and !failure.work_complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    try std.testing.expectEqual(.not_attempted, failure.outputs.stdout.final_flush);
    try std.testing.expectEqual(@as(usize, 0), failure.secondaryItems().len);
}

fn expectComplete(failure: *const diagnostics.Failure) !void {
    try std.testing.expect(failure.primary() == null);
    try std.testing.expect(failure.complete and failure.work_complete);
    try std.testing.expectEqual(.succeeded, failure.outputs.stdout.final_flush);
}
