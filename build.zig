//! Builds the `z-fasta` executable, the `z-fasta` package (`src/root.zig`), and the tests.
//! Command tests run in a build-owned copy of the declared fixtures, never in `tests/data/`.

const std = @import("std");
const package = @import("build.zig.zon");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const test_filter = b.option([]const u8, "test-filter", "Run tests whose names contain this text");
    const test_filters: []const []const u8 = if (test_filter) |filter| &.{filter} else &.{};
    const strip = b.option(bool, "strip", "Strip debug info") orelse (optimize != .Debug);

    const diagnostics_module = b.createModule(.{
        .root_source_file = b.path("src/diagnostics.zig"),
        .target = target,
        .optimize = optimize,
    });
    const library_module = b.addModule("z-fasta", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
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

    const build_options = b.addOptions();
    build_options.addOption([]const u8, "version", package.version);
    const build_options_module = build_options.createModule();

    exe.root_module.addImport("diagnostics", diagnostics_module);
    exe.root_module.addImport("build_options", build_options_module);
    if (optimize == .ReleaseSmall) {
        const unsupported = b.addFail("ReleaseSmall is unsupported for the z-fasta executable; use -Doptimize=ReleaseFast (strips by default)");
        exe.step.dependOn(&unsupported.step);
    }
    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run z-fasta");
    run_step.dependOn(&run_cmd.step);

    const main_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    main_module.addImport("diagnostics", diagnostics_module);
    main_module.addImport("build_options", build_options_module);
    const test_options = b.addOptions();
    test_options.addOptionPath("binary_path", exe.getEmittedBin());
    test_options.addOption([]const u8, "build_root", b.pathFromRoot("."));
    test_options.addOption([]const u8, "version", package.version);
    const test_options_module = test_options.createModule();

    const test_index_module = b.createModule(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "main", .module = main_module },
            .{ .name = "test_options", .module = test_options_module },
        },
    });
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

    const run_test_getter_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/getter.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
        }),
    }));

    const run_test_main_unit = b.addRunArtifact(b.addTest(.{ .filters = test_filters, .root_module = main_module }));

    const run_test_library = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_library.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "z-fasta", .module = library_module },
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));

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

    const run_test_stats_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/stats.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
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

    const run_test_fasta_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/fasta.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "diagnostics", .module = diagnostics_module }},
        }),
    }));

    const run_test_index_format_unit = b.addRunArtifact(b.addTest(.{
        .filters = test_filters,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/index_format.zig"),
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
                .{ .name = "test_options", .module = test_options_module },
            },
        }),
    }));
    const run_consumer_diagnostics = b.addRunArtifact(b.addExecutable(.{
        .name = "consumer-diagnostics",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/consumer_diagnostics.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "diagnostics", .module = diagnostics_module },
                .{ .name = "reporting", .module = reporting_module },
            },
        }),
    }));
    const run_consumer_commands = b.addRunArtifact(b.addExecutable(.{
        .name = "consumer-commands",
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/consumer_commands.zig"),
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
    diagnostics_step.dependOn(&run_consumer_commands.step);
    diagnostics_step.dependOn(&run_test_diagnostics.step);
    diagnostics_step.dependOn(&run_consumer_diagnostics.step);

    const library_bench = b.addExecutable(.{
        .name = "z-fasta-library",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/get/library.zig"),
            .target = target,
            .optimize = optimize,
            .strip = strip,
            .imports = &.{.{ .name = "z-fasta", .module = library_module }},
        }),
    });
    const library_bench_step = b.step("bench-library", "Build the package lane of the GET benchmark");
    library_bench_step.dependOn(&b.addInstallArtifact(library_bench, .{}).step);

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
        run_test_index,        run_test_get,         run_test_main_unit,   run_test_indexer_unit,
        run_test_stats,        run_test_stats_unit,  run_test_validate,    run_test_validator_unit,
        run_consumer_commands, run_test_diagnostics, run_test_getter_unit,
    }) |run| run.setCwd(fixture_root);

    const index_step = b.step("test-index", "Run index tests");
    index_step.dependOn(&run_test_index.step);
    index_step.dependOn(&run_test_indexer_unit.step);
    index_step.dependOn(&run_test_index_format_unit.step);
    const get_step = b.step("test-get", "Run GET tests");
    get_step.dependOn(&run_test_get.step);
    get_step.dependOn(&run_test_getter_unit.step);
    const stats_step = b.step("test-stats", "Run stats tests");
    stats_step.dependOn(&run_test_stats.step);
    stats_step.dependOn(&run_test_stats_unit.step);
    const validate_step = b.step("test-validate", "Run validator tests");
    validate_step.dependOn(&run_test_validate.step);
    validate_step.dependOn(&run_test_validator_unit.step);
    const cli_step = b.step("test-cli", "Run CLI dispatch tests");
    cli_step.dependOn(&run_test_main_unit.step);
    cli_step.dependOn(&run_test_diagnostics.step);
    const library_step = b.step("test-library", "Run public package tests");
    library_step.dependOn(&run_test_library.step);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(diagnostics_step);
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_getter_unit.step);
    test_step.dependOn(&run_test_main_unit.step);
    test_step.dependOn(&run_test_indexer_unit.step);
    test_step.dependOn(&run_test_stats.step);
    test_step.dependOn(&run_test_stats_unit.step);
    test_step.dependOn(&run_test_validate.step);
    test_step.dependOn(&run_test_validator_unit.step);
    test_step.dependOn(&run_test_fasta_unit.step);
    test_step.dependOn(&run_test_index_format_unit.step);
    test_step.dependOn(&run_test_library.step);
}
