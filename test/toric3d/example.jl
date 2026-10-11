using Test

const TORIC3D_EXAMPLE = abspath(joinpath(
    @__DIR__, "..", "..", "examples", "compare_toric3d_feedback.jl"))
const TORIC3D_PROJECT = abspath(joinpath(@__DIR__, "..", ".."))

function _toric3d_example_process(arguments...)
    command = `$(Base.julia_cmd()) --project=$(TORIC3D_PROJECT) $(TORIC3D_EXAMPLE) $(collect(arguments))`
    mktemp() do path, io
        process = run(pipeline(ignorestatus(command); stdout=io, stderr=io))
        seekstart(io)
        return (success(process), read(io, String))
    end
end

@testset "3D feedback example CLI" begin
    @test isfile(TORIC3D_EXAMPLE)

    ok, output = _toric3d_example_process("--help")
    @test ok
    @test contains(output, "--pilot-shots")
    @test contains(output, "--output-dir")

    for arguments in (("--sizes", "3,2"), ("--shots", "0"),
                      ("--error-rates", "0,0"), ("--pilot-shots", "0"))
        ok, output = _toric3d_example_process(arguments...)
        @test !ok
        @test contains(output, "error:")
        @test contains(output, "Usage:")
    end

    mktempdir() do directory
        ok, output = _toric3d_example_process(
            "--sizes", "2", "--error-rates", "0", "--shots", "1",
            "--pilot-shots", "1", "--seed", "1234",
            "--output-dir", directory)
        @test ok
        @test contains(output, "Warmup including compilation:")
        @test contains(output, "Pilot runtime:")
        @test contains(output, "Rough full-scan estimate:")
        @test contains(output, "paired shots in")
        @test contains(output, "Saved comparison:")
        @test isfile(joinpath(directory, "toric3d-feedback-raw.csv"))
        @test isfile(joinpath(directory, "toric3d-feedback-metadata.toml"))
        for extension in ("svg", "pdf", "png")
            @test isfile(joinpath(directory, "toric3d-feedback.$extension"))
        end
    end
end
