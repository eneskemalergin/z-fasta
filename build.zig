const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "z-fasta",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    b.installArtifact(exe);

    // Run command: zig build run -- <args>
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run z-fasta");
    run_step.dependOn(&run_cmd.step);

    // Test command: zig build test
    const main_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
    });

    // Index tests
    const test_index = b.addTest(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
    });
    test_index.root_module.addImport("main", main_module);
    const run_test_index = b.addRunArtifact(test_index);

    // Get tests
    const test_get = b.addTest(.{
        .root_source_file = b.path("tests/test_get.zig"),
        .target = target,
        .optimize = optimize,
    });
    test_get.root_module.addImport("main", main_module);
    const run_test_get = b.addRunArtifact(test_get);

    // Stats tests
    const test_stats = b.addTest(.{
        .root_source_file = b.path("tests/test_stats.zig"),
        .target = target,
        .optimize = optimize,
    });
    test_stats.root_module.addImport("main", main_module);
    const run_test_stats = b.addRunArtifact(test_stats);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_stats.step);
}
