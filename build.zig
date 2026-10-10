//! Builds the z-fasta executable and its unit, integration, and CLI test suites.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const test_filter = b.option([]const u8, "test-filter", "Run tests whose names contain this text");
    const test_filters: []const []const u8 = if (test_filter) |filter| &.{filter} else &.{};
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

    const diagnostics_module = b.createModule(.{
        .root_source_file = b.path("src/diagnostics.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe = b.addExecutable(.{
        .name = "z-fasta",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .strip = strip,
        }),
    });

    exe.root_module.addImport("diagnostics", diagnostics_module);
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

    main_module.addImport("diagnostics", diagnostics_module);
    const test_options = b.addOptions();
    test_options.addOptionPath("binary_path", exe.getEmittedBin());
    test_options.addOption([]const u8, "build_root", b.pathFromRoot("."));
    const test_options_module = test_options.createModule();
    main_module.addImport("test_options", test_options_module);

    const test_index_module = b.createModule(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "main", .module = main_module },
            .{ .name = "test_options", .module = test_options_module },
        },
    });
    // Timestamp tests call POSIX time functions through std.c.
    test_index_module.link_libc = true;
    const run_test_index = b.addRunArtifact(b.addTest(.{ .filters = test_filters, .root_module = test_index_module }));

    const run_test_get = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_get.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "main", .module = main_module },
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));

    const run_test_main_unit = b.addRunArtifact(b.addTest(.{ .filters = test_filters, .root_module = main_module }));

    const run_test_indexer_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/indexer.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
        }),
    }));

    const run_test_stats = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_stats.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "main", .module = main_module },
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));

    const run_test_validator_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/validator.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
        }),
    }));

    const run_test_validate = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_validate.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "main", .module = main_module },
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));

    const reporting_module = b.createModule(.{
        .root_source_file = b.path("src/reporting.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
    });
    const run_test_reporting = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_reporting.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "diagnostics", .module = diagnostics_module },
                .{ .name = "reporting", .module = reporting_module },
            },
        }),
    }));
    const reporting_step = b.step("test-reporting", "Run diagnostic renderer tests");
    reporting_step.dependOn(&run_test_reporting.step);
    const run_test_diagnostics = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_diagnostics.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "diagnostics", .module = diagnostics_module },
                .{ .name = "reporting", .module = reporting_module },
            },
        }),
    }));
    const run_diagnostics_caller = b.addRunArtifact(b.addExecutable(.{
        .name = "diagnostics-caller",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/diagnostics_caller.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "diagnostics", .module = diagnostics_module },
                .{ .name = "reporting", .module = reporting_module },
            },
        }),
    }));
    const run_commands_caller = b.addRunArtifact(b.addExecutable(.{
        .name = "diagnostics-commands-caller",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/diagnostics_commands_caller.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "main", .module = main_module },
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));
    const diagnostics_step = b.step("test-diagnostics", "Run diagnostics tests and separate callers");
    diagnostics_step.dependOn(reporting_step);
    diagnostics_step.dependOn(&run_commands_caller.step);
    diagnostics_step.dependOn(&run_test_diagnostics.step);
    diagnostics_step.dependOn(&run_diagnostics_caller.step);

    const prepare_fixtures = b.addRunArtifact(b.addExecutable(.{
        .name = "prepare-test-fixtures",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/prepare_fixtures.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    const fixture_root = prepare_fixtures.addOutputDirectoryArg("test-fixtures");
    // Only declared fixtures enter the test directory, excluding local sidecars.
    for ([_][]const u8{
        "tests/data/edge_cases.fasta",
        "tests/data/empty.fasta",
        "tests/data/gates/duplicates.fasta",
        "tests/data/gates/invalid_utf8_header.fasta",
        "tests/data/gates/long_header.fasta",
        "tests/data/gates/long_header_reject.fasta",
        "tests/data/gates/seq1_large_offset.fasta",
        "tests/data/gates/seq1_large_offset.fasta.fai",
        "tests/data/gates/seq1_zero_geometry.fasta",
        "tests/data/gates/seq1_zero_geometry.fasta.fai",
        "tests/data/mixed_widths.fasta",
        "tests/data/not_fasta.txt",
        "tests/data/proteome.fasta",
        "tests/data/simple.fasta",
        "tests/data/single.fasta",
        "tests/data/validator_indexer_agreement.fasta",
    }) |path| {
        prepare_fixtures.addArg(path);
        prepare_fixtures.addFileArg(b.path(path));
    }
    for ([_]*std.Build.Step.Run{
        run_test_index, run_test_get,      run_test_main_unit,      run_test_indexer_unit,
        run_test_stats, run_test_validate, run_test_validator_unit, run_commands_caller,
    }) |run| run.setCwd(fixture_root);

    const index_step = b.step("test-index", "Run index tests");
    index_step.dependOn(&run_test_index.step);
    index_step.dependOn(&run_test_indexer_unit.step);
    const get_step = b.step("test-get", "Run GET tests");
    get_step.dependOn(&run_test_get.step);
    const stats_step = b.step("test-stats", "Run stats tests");
    stats_step.dependOn(&run_test_stats.step);
    const validate_step = b.step("test-validate", "Run validator tests");
    validate_step.dependOn(&run_test_validate.step);
    validate_step.dependOn(&run_test_validator_unit.step);
    const cli_step = b.step("test-cli", "Run CLI dispatch tests");
    cli_step.dependOn(&run_test_main_unit.step);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(diagnostics_step);
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_main_unit.step);
    test_step.dependOn(&run_test_indexer_unit.step);
    test_step.dependOn(&run_test_stats.step);
    test_step.dependOn(&run_test_validate.step);
    test_step.dependOn(&run_test_validator_unit.step);
}
