using Test
using TopoNoise

const TORIC3D_DIAGONAL_FEEDBACK_EXAMPLE = abspath(joinpath(
    @__DIR__, "..", "..", "examples", "compare_toric3d_diagonal_feedback.jl"))
include(TORIC3D_DIAGONAL_FEEDBACK_EXAMPLE)

function _diagonal_feedback_main(arguments)
    output = IOBuffer()
    errors = IOBuffer()
    status = Toric3DDiagonalFeedbackExample.main(
        arguments; io=output, error_io=errors)
    return status, String(take!(output)), String(take!(errors))
end

@testset "Diagonal three-feedback example CLI" begin
    status, output, errors = _diagonal_feedback_main(["--help"])
    @test status == 0
    @test isempty(errors)
    for option in ("--sizes", "--rates", "--channels", "--shots",
            "--batch-size", "--seed", "--x-routing", "--output",
            "--resume", "--timing-shots", "--pilot-only")
        @test contains(output, option)
    end

    x_defaults = Toric3DDiagonalFeedbackExample._parse_arguments(String[])
    @test x_defaults.channels == [:x_only]
    @test x_defaults.sizes == [3, 5, 7, 9]
    @test x_defaults.shots == 200

    z_defaults = Toric3DDiagonalFeedbackExample._parse_arguments(
        ["--channels", "z_only"])
    @test z_defaults.channels == [:z_only]
    @test z_defaults.sizes == [5, 9, 13, 17]

    explicit = Toric3DDiagonalFeedbackExample._parse_arguments([
        "--sizes", "1x2x2,2", "--rates", "0,0.01",
        "--channels", "x_only,z_only", "--shots", "3",
        "--batch-size", "2", "--seed", "7", "--x-routing", "full_path",
        "--output", "out", "--timing-shots", "1", "--resume",
        "--pilot-only"])
    @test explicit.sizes == [(1, 2, 2), (2, 2, 2)]
    @test explicit.error_rates == [0.0, 0.01]
    @test explicit.channels == [:x_only, :z_only]
    @test explicit.resume && explicit.pilot_only
    @test explicit.x_routing === :full_path

    invalid = [
        ["--unknown", "1"],
        ["--sizes"],
        ["--sizes", "2,1"],
        ["--rates", "0.01,0"],
        ["--channels", "bad"],
        ["--channels", "x_only,z_only"],
        ["--channels", "x_only,z_only", "--sizes", "2"],
        ["--shots", "0"],
        ["--batch-size", "0"],
        ["--timing-shots", "0"],
        ["--seed", "-1"],
        ["--x-routing", "bad"],
        ["--output", ""]]
    for arguments in invalid
        status, _, errors = _diagonal_feedback_main(arguments)
        @test status == 1
        @test !isempty(errors)
        @test contains(errors, "Usage:")
    end

    mktempdir() do directory
        arguments = ["--sizes", "1x2x2", "--rates", "0",
            "--channels", "x_only", "--shots", "2", "--batch-size", "1",
            "--timing-shots", "1", "--seed", "1234", "--x-routing",
            "one_step", "--output", directory]
        status, output, errors = _diagonal_feedback_main(arguments)
        @test status == 0
        @test isempty(errors)
        @test contains(output, "Warmup")
        @test contains(output, "Timing pilot")
        @test contains(output, "Rough scan estimate")
        @test contains(output, "Saved")

        expected = ["timing-pilot.toml",
            "toric3d-diagonal-feedback-raw.csv",
            "toric3d-diagonal-feedback-metadata.toml",
            "toric3d-diagonal-feedback.svg", "toric3d-diagonal-feedback.pdf",
            "toric3d-diagonal-feedback.png",
            "toric3d-diagonal-feedback-checkpoint.toml"]
        @test all(name -> isfile(joinpath(directory, name)), expected)

        timing = TopoNoise.TOML.parsefile(joinpath(directory, "timing-pilot.toml"))
        @test timing["timing_shots"] == 1
        @test timing["production_shots_per_point"] == 2
        @test !timing["timing_shots_included_in_production"]
        @test timing["decoder_construction_seconds"] >= 0
        for protocol in ("none", "slice_local", "full_prefix")
            @test timing["protocol_seconds"][protocol] >= 0
        end

        rows = readlines(joinpath(directory,
            "toric3d-diagonal-feedback-raw.csv"))
        @test length(rows) == 4
        header = split(first(rows), ',')
        first_joint = findfirst(==("joint_n0_s0_f0"), header)
        last_joint = findfirst(==("joint_n1_s1_f1"), header)
        @test all(sum(parse.(Int, split(row, ',')[first_joint:last_joint])) == 2
            for row in rows[2:end])
    end

    mktempdir() do directory
        arguments = ["--sizes", "1x2x2", "--rates", "0",
            "--channels", "x_only", "--shots", "2", "--batch-size", "1",
            "--timing-shots", "1", "--seed", "1234", "--output", directory,
            "--pilot-only"]
        status, output, errors = _diagonal_feedback_main(arguments)
        @test status == 0
        @test isempty(errors)
        @test contains(output, "Pilot-only")
        @test isfile(joinpath(directory, "timing-pilot.toml"))
        @test !isfile(joinpath(directory, "toric3d-diagonal-feedback-raw.csv"))
        @test !isfile(joinpath(directory,
            "toric3d-diagonal-feedback-checkpoint.toml"))
    end
end
