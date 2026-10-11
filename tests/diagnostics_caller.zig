//! Separate program that uses `diagnostics` directly, without the CLI, to check that a caller can
//! handle failures and release partial allocations. Faults are injected, not real disk errors.

const std = @import("std");
const diagnostics = @import("diagnostics");
const reporting = @import("reporting");

pub fn exercise(allocator: std.mem.Allocator) !void {
    var failure: diagnostics.Failure = undefined;
    const caller_path = [_]u8{ 's', 'a', 'm', 'p', 'l', 'e', '.', 'f', 'a' };
    var failing = std.testing.FailingAllocator.init(allocator, .{ .fail_index = 1 });
    if (acquire(failing.allocator(), &failure, &caller_path, false)) |owned| {
        owned.deinit(failing.allocator());
        return error.ExpectedAllocationFailure;
    } else |err| {
        try require(err == error.OutOfMemory);
        try require(failure.primary().?.code() == .@"resource.out_of_memory");
        const facts = failure.primary().?.details.@"resource.out_of_memory";
        try require(facts.allocation_role == .records and facts.requested_bytes.? == 128);
    }
    try require(failing.allocations == 1 and failing.deallocations == 1);
    try require(failing.allocated_bytes == failing.freed_bytes);
    const source = failure.sources.get(.zfi).?;
    try require(source.label.borrowed.ptr == &caller_path and source.suffix == .zfi);

    if (acquire(allocator, &failure, &caller_path, true)) |owned| {
        owned.deinit(allocator);
        return error.ExpectedParseFailure;
    } else |err| {
        try require(err == error.InvalidIndex);
        try require(failure.primary().?.code() == .@"index.fai.invalid_integer");
        const retained = failure.primary().?.details.@"index.fai.invalid_integer".value.?;
        try require(std.mem.eql(u8, failure.preview(retained), "bad\xff!"));
        try require(failure.secondaryItems().len == 1);
        try require(failure.secondaryItems()[0].cause.?.zig == error.WriteFailed);
        try require(failure.outputs.file.state == .staged);
    }
    try require(failure.sources.get(.fasta).?.label.borrowed.ptr == &caller_path);

    var report_buffer: [reporting.JSON_BYTES_LIMIT]u8 = undefined;
    var report_writer = std.Io.Writer.fixed(&report_buffer);
    try reporting.writeJson(&report_writer, &failure, .{});
    try require(std.mem.indexOf(u8, report_writer.buffered(), "YmFk/yE=") != null);
    report_writer.end = 0;
    try reporting.writeText(&report_writer, &failure);
    try require(std.mem.indexOf(u8, report_writer.buffered(), "bad\\xff!") != null);

    const owned = try acquire(allocator, &failure, &caller_path, false);
    defer owned.deinit(allocator);
    try require(failure.primary() == null and failure.secondary_len == 0);
    try require(failure.preview_len == 0 and !failure.secondary_truncated);
    try require(failure.complete and failure.work_complete);
    try require(failure.outputs.file.state == .not_started);
    try require(owned.name.len == 64 and owned.records.len == 128);
}

pub fn main() !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    const result = exercise(debug_allocator.allocator());
    const cleanup = debug_allocator.deinit();
    if (cleanup == .leak) return error.CallerMemoryLeak;
    try result;
}

const Owned = struct {
    name: []u8,
    records: []u8,

    fn deinit(self: Owned, allocator: std.mem.Allocator) void {
        allocator.free(self.records);
        allocator.free(self.name);
    }
};

fn acquire(
    allocator: std.mem.Allocator,
    failure: *diagnostics.Failure,
    caller_path: []const u8,
    reject_after_progress: bool,
) diagnostics.Error!Owned {
    failure.reset(.get);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = caller_path }, .kind = .path });
    failure.sources.set(.zfi, .{ .label = .{ .borrowed = caller_path }, .kind = .path, .suffix = .zfi });
    failure.sources.set(.fai, .{ .label = .{ .borrowed = caller_path }, .kind = .path, .suffix = .fai });
    const name = allocator.alloc(u8, 64) catch {
        @branchHint(.cold);
        return failure.fail(.{
            .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .names, .requested_bytes = 64 } },
            .stage = .allocate,
            .cause = .{ .zig = error.OutOfMemory },
        });
    };
    errdefer allocator.free(name);
    const records = allocator.alloc(u8, 128) catch {
        @branchHint(.cold);
        return failure.fail(.{
            .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .records, .requested_bytes = 128 } },
            .stage = .allocate,
            .cause = .{ .zig = error.OutOfMemory },
        });
    };
    errdefer allocator.free(records);

    if (reject_after_progress) {
        @branchHint(.cold);
        @memcpy(name[0..5], "bad\xff!");
        const preview = failure.retain(name[0..5], diagnostics.PREVIEW_BYTES);
        var buffer: [1]u8 = undefined;
        var writer = std.Io.Writer.fixed(&buffer);
        writer.writeAll("staged") catch |err| {
            failure.deferFailure(.{
                .details = .{ .@"io.write_failed" = {} },
                .stage = .write,
                .cause = .{ .zig = err },
            });
        };
        failure.outputs.file.state = .staged;
        return failure.fail(.{
            .details = .{ .@"index.fai.invalid_integer" = .{ .field = .length, .value = preview } },
            .stage = .parse,
            .at = .{ .source = .fai, .location = .{ .offset = 0, .line = 1, .column = 2 } },
        });
    }
    failure.succeed();
    return .{ .name = name, .records = records };
}

fn require(condition: bool) error{CallerCheckFailed}!void {
    if (!condition) return error.CallerCheckFailed;
}
