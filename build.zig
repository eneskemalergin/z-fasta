const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "z-fasta",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
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

    // Shared main module for test imports
    const main_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Index tests
    const test_index_module = b.createModule(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "main", .module = main_module }},
    });
    test_index_module.link_libc = true;
    const test_index = b.addTest(.{
        .root_module = test_index_module,
    });
    const run_test_index = b.addRunArtifact(test_index);

    // Get tests
    const test_get = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_get.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    });
    const run_test_get = b.addRunArtifact(test_get);

    // Stats tests
    const test_stats = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_stats.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    });
    const run_test_stats = b.addRunArtifact(test_stats);

    const test_complement = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/complement.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_test_complement = b.addRunArtifact(test_complement);

    const test_bed_parser = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/bed_parser.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_test_bed_parser = b.addRunArtifact(test_bed_parser);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_stats.step);
    test_step.dependOn(&run_test_complement.step);
    test_step.dependOn(&run_test_bed_parser.step);
}
