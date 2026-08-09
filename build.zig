//! Builds the z-fasta executable and its unit, integration, and CLI test suites.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    // ReleaseSmall regresses the text hot path; ReleaseFast is the supported ship mode.
    if (optimize == .ReleaseSmall) {
        std.debug.print(
            "error: ReleaseSmall is unsupported; use -Doptimize=ReleaseFast (strips by default)\n",
            .{},
        );
        std.process.exit(1);
    }
    // ReleaseFast does not enable stripping automatically.
    const strip = b.option(bool, "strip", "Strip debug info") orelse (optimize != .Debug);

    const exe = b.addExecutable(.{
        .name = "z-fasta",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .strip = strip,
        }),
    });

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run z-fasta");
    run_step.dependOn(&run_cmd.step);

    // Test modules stay unstripped; ship strip applies to the installed exe only.
    const main_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    const test_index_module = b.createModule(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "main", .module = main_module }},
    });
    // Timestamp tests call POSIX time functions through std.c.
    test_index_module.link_libc = true;
    const run_test_index = b.addRunArtifact(b.addTest(.{ .root_module = test_index_module }));
    // CLI-bearing suites spawn the installed path, so keep it current.
    run_test_index.step.dependOn(b.getInstallStep());

    const run_test_get = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_get.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    run_test_get.step.dependOn(b.getInstallStep());

    const run_test_main_unit = b.addRunArtifact(b.addTest(.{ .root_module = main_module }));
    run_test_main_unit.step.dependOn(b.getInstallStep());

    const run_test_indexer_unit = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/indexer.zig"),
            .target = target,
            .optimize = optimize,
        }),
    }));

    const run_test_stats = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_stats.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    run_test_stats.step.dependOn(b.getInstallStep());

    // Generate ignored sidecars so tests also pass on a clean checkout.
    const gen_simple_index = b.addRunArtifact(exe);
    gen_simple_index.addArgs(&.{ "index", "tests/data/simple.fasta" });
    run_test_index.step.dependOn(&gen_simple_index.step);
    run_test_get.step.dependOn(&gen_simple_index.step);
    run_test_main_unit.step.dependOn(&gen_simple_index.step);
    run_test_stats.step.dependOn(&gen_simple_index.step);

    const gen_proteome_index = b.addRunArtifact(exe);
    gen_proteome_index.addArgs(&.{ "index", "tests/data/proteome.fasta" });
    run_test_get.step.dependOn(&gen_proteome_index.step);
    run_test_stats.step.dependOn(&gen_proteome_index.step);

    const gen_edge_cases_index = b.addRunArtifact(exe);
    gen_edge_cases_index.addArgs(&.{ "index", "tests/data/edge_cases.fasta" });
    run_test_get.step.dependOn(&gen_edge_cases_index.step);

    const run_test_validate = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_validate.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    run_test_validate.step.dependOn(b.getInstallStep());

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_main_unit.step);
    test_step.dependOn(&run_test_indexer_unit.step);
    test_step.dependOn(&run_test_stats.step);
    test_step.dependOn(&run_test_validate.step);
}
