//! Separate executable calling the four core operations through existing exports.
//! No CLI dispatcher runs. The caller handles failures, then calls each command
//! again with valid input, discarding only stdout through the test backend.

const std = @import("std");
const commands = @import("main");

pub fn main(init: std.process.Init.Minimal) !void {
    const allocator = std.heap.page_allocator;
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    backend = threaded.io();
    var vtable = backend.vtable.*;
    vtable.operate = discardStdout;
    const io = std.Io{ .userdata = backend.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    const missing = "tests/data/diagnostics-caller-does-not-exist.fa";
    const valid = "tests/data/simple.fasta";
    try std.testing.expectError(error.Io, commands.indexer.runIndex(io, init.environ, missing, .{ .emit_fai = true }, &failure));
    try expectMissing(&failure, .index);
    const indexed = try commands.indexer.runIndex(io, init.environ, valid, .{ .emit_fai = true }, &failure);
    try std.testing.expectEqual(@as(u64, 2), indexed.record_count);
    try expectComplete(&failure);

    const get_options = commands.getter.Options{ .source = .{ .positional = &.{"seq1:1-4"} } };
    try std.testing.expectError(error.Io, commands.getter.runGetWithOptions(allocator, io, missing, get_options, &failure));
    try expectMissing(&failure, .get);
    const extracted = try commands.getter.runGetWithOptions(allocator, io, valid, get_options, &failure);
    try std.testing.expectEqual(@as(u64, 4), extracted.total_bases);
    try expectComplete(&failure);

    try std.testing.expectError(error.Io, commands.stats.runStats(allocator, io, missing, &failure));
    try expectMissing(&failure, .stats);
    try commands.stats.runStats(allocator, io, valid, &failure);
    try expectComplete(&failure);

    try std.testing.expectError(error.Io, commands.validator.runValidate(allocator, io, missing, .{}, &failure));
    try expectMissing(&failure, .validate);
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

fn expectMissing(failure: *const diagnostics.Failure, operation: diagnostics.Operation) !void {
    try std.testing.expectEqual(operation, failure.operation);
    try std.testing.expectEqual(diagnostics.Code.@"io.file_not_found", failure.primary().?.code());
    try std.testing.expectEqual(error.FileNotFound, failure.primary().?.cause.?.zig);
    try std.testing.expect(!failure.complete and !failure.work_complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
}

fn expectComplete(failure: *const diagnostics.Failure) !void {
    try std.testing.expect(failure.primary() == null);
    try std.testing.expect(failure.complete and failure.work_complete);
    try std.testing.expectEqual(.succeeded, failure.outputs.stdout.final_flush);
}
