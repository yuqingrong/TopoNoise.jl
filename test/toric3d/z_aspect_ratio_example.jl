using Test
using TopoNoise

const TORIC3D_Z_ASPECT_RATIO_EXAMPLE = abspath(joinpath(
    @__DIR__, "..", "..", "examples", "compare_toric3d_z_aspect_ratio.jl"))
include(TORIC3D_Z_ASPECT_RATIO_EXAMPLE)

function _z_aspect_main(arguments)
    output = IOBuffer()
    errors = IOBuffer()
    status = Toric3DZAspectRatioExample.main(arguments; io=output, error_io=errors)
    return status, String(take!(output)), String(take!(errors))
end

@testset "H+6 Z aspect-ratio example CLI" begin
    status, output, errors = _z_aspect_main(["--help"])
    @test status == 0
    @test isempty(errors)
    @test contains(output, "--step-multipliers")
    @test contains(output, "final-only")
    @test contains(output, "no intermediate feedback")

    defaults = Toric3DZAspectRatioExample._parse_arguments(String[])
    @test defaults.sizes == [3,5,7]
    @test defaults.step_multipliers == [1,2]
    @test defaults.error_rates == [.0045,.005,.0055,.006,.0065]
    @test defaults.shots == 1000
    @test defaults.seed == 1234
    @test defaults.batch_size == 256
    @test defaults.pilot_shots == 20
    @test defaults.output_dir ==
        "results/toric3d-h6-z-aspect-ratio-2026-10-15/production"
    @test defaults.dimensions == [
        (3,3,3),(6,3,3),(5,5,5),(10,5,5),(7,7,7),(14,7,7)]

    parsed = Toric3DZAspectRatioExample._parse_arguments([
        "--sizes", "2,4", "--step-multipliers", "3,1",
        "--error-rates", "0,.01", "--shots", "2", "--seed", "9",
        "--batch-size", "1", "--pilot-shots", "1", "--output-dir", "out"])
    @test parsed.dimensions == [(6,2,2),(2,2,2),(12,4,4),(4,4,4)]
    @test_throws ArgumentError Toric3DZAspectRatioExample._aspect_dimensions([3.0],[1])
    @test_throws ArgumentError Toric3DZAspectRatioExample._aspect_dimensions([3],[1.0])
    @test Toric3DZAspectRatioExample._timing_target(
        [(8,2,2),(2,5,5)]) == (2,5,5)

    invalid = [
        ["--unknown", "1"],
        ["--sizes"],
        ["--sizes", ""],
        ["--sizes", "0"],
        ["--sizes", "1"],
        ["--sizes", "3,3"],
        ["--step-multipliers", "0"],
        ["--step-multipliers", "2,2"],
        ["--error-rates", "-0.01"],
        ["--error-rates", "0.5"],
        ["--error-rates", "0.01,0.01"],
        ["--shots", "0"],
        ["--seed", "-1"],
        ["--batch-size", "0"],
        ["--pilot-shots", "0"],
        ["--output-dir", ""]]
    for arguments in invalid
        status, _, errors = _z_aspect_main(arguments)
        @test status == 1
        @test contains(errors, "Usage:")
    end

    mktempdir() do directory
        arguments = ["--sizes", "2", "--step-multipliers", "1",
            "--error-rates", "0", "--shots", "1", "--batch-size", "1",
            "--pilot-shots", "1", "--seed", "1234", "--output-dir", directory]
        status, output, errors = _z_aspect_main(arguments)
        @test status == 0
        @test isempty(errors)
        @test contains(output, "Timing pilot")
        @test contains(output, "Saved")

        expected = [
            "timing-pilot.toml",
            "toric3d-z-threshold-raw.csv",
            "toric3d-z-threshold-metadata.toml",
            "toric3d-z-aspect-ratio-study.toml"]
        @test all(name -> isfile(joinpath(directory,name)), expected)

        timing = TopoNoise.TOML.parsefile(joinpath(directory,"timing-pilot.toml"))
        @test timing["pilot_dimensions_steps_Ly_Lz"] == [2,2,2]
        @test timing["pilot_shots"] == 1
        @test timing["production_shots_per_point"] == 1
        @test !timing["timing_shots_included_in_production"]

        study = TopoNoise.TOML.parsefile(joinpath(
            directory,"toric3d-z-aspect-ratio-study.toml"))
        @test study["circuit"] == "plaquette-transfer H+6-CNOT"
        @test study["decoding"] == "final-only"
        @test study["intermediate_feedback"] == "none"
        @test study["comparison"] == "T=L versus T=2L aspect-ratio diagnostic"
        @test study["dimensions_steps_Ly_Lz"] == [[2,2,2]]

        rows = readlines(joinpath(directory,"toric3d-z-threshold-raw.csv"))
        @test length(rows) == 2
        @test split(rows[2],',')[6] == "0"
    end
end
