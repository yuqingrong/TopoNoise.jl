using Test
using TopoNoise

const TORIC3D_X_FEEDBACK_EXAMPLE = abspath(joinpath(
    @__DIR__, "..", "..", "examples", "scan_toric3d_x_feedback_threshold.jl"))
include(TORIC3D_X_FEEDBACK_EXAMPLE)

function _x_feedback_main(arguments)
    output = IOBuffer()
    errors = IOBuffer()
    status = Toric3DXFeedbackExample.main(arguments; io=output, error_io=errors)
    return status, String(take!(output)), String(take!(errors))
end

@testset "Diagonal X-feedback example CLI" begin
    status, output, errors = _x_feedback_main(["--help"])
    @test status == 0
    @test isempty(errors)
    @test contains(output, "--resume")
    @test contains(output, "--batch-size")
    @test contains(output, "--routing one_step|full_path")

    parsed = Toric3DXFeedbackExample._parse_arguments([
        "--resume", "--sizes", "2", "--error-rates", "0",
        "--shots", "1", "--batch-size", "1"])
    @test parsed.resume
    @test parsed.batch_size == 1

    invalid = [
        ["--unknown", "1"],
        ["--sizes"],
        ["--sizes", "3,2"],
        ["--error-rates", "0.01,0"],
        ["--shots", "0"],
        ["--batch-size", "0"],
        ["--pilot-shots", "0"],
        ["--routing", "bad"],
        ["--output-dir", ""]]
    for arguments in invalid
        status, _, errors = _x_feedback_main(arguments)
        @test status == 1
        @test !isempty(errors)
        @test contains(errors, "Usage:")
    end

    mktempdir() do directory
        arguments = ["--sizes", "2", "--error-rates", "0", "--shots", "1",
            "--batch-size", "1", "--pilot-shots", "1", "--seed", "1234",
            "--routing", "one_step", "--output-dir", directory]
        status, output, errors = _x_feedback_main(arguments)
        @test status == 0
        @test isempty(errors)
        @test contains(output, "Warmup")
        @test contains(output, "Timing pilot")
        @test contains(output, "Rough scan estimate")
        @test contains(output, "Saved")

        expected = [
            "timing-pilot.toml",
            "toric3d-x-feedback-raw.csv",
            "toric3d-x-feedback-metadata.toml",
            "toric3d-x-feedback.svg",
            "toric3d-x-feedback.pdf",
            "toric3d-x-feedback.png",
            "toric3d-x-feedback-checkpoint.toml"]
        @test all(name -> isfile(joinpath(directory, name)), expected)

        timing = TopoNoise.TOML.parsefile(joinpath(directory, "timing-pilot.toml"))
        @test timing["warmup_shots"] == 2
        @test timing["pilot_shots"] == 1
        @test timing["production_shots_per_point"] == 1
        @test !timing["timing_shots_included_in_production"]
        @test timing["setup_seconds"] >= 0
        @test timing["first_decode_seconds"] >= 0
        @test timing["pilot_seconds"] >= 0

        rows = readlines(joinpath(directory, "toric3d-x-feedback-raw.csv"))
        @test length(rows) == 4

        status, output, errors = _x_feedback_main(vcat(arguments, ["--resume"]))
        @test status == 0
        @test isempty(errors)
        @test contains(output, "Resuming")
    end
end
