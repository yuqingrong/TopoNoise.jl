using Test

const _OPEN_CODE_CAPACITY_EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "examples", "scan_open_code_capacity.jl")

@testset "open code-capacity example" begin
    @test isfile(_OPEN_CODE_CAPACITY_EXAMPLE_PATH)

    if isfile(_OPEN_CODE_CAPACITY_EXAMPLE_PATH)
        include(_OPEN_CODE_CAPACITY_EXAMPLE_PATH)

        mktempdir() do directory
            output = IOBuffer()
            status = ScanOpenCodeCapacityExample.main([
                "--p-min", "0.05", "--p-max", "0.10", "--p-step", "0.05",
                "--sizes", "3,4", "--shots", "12", "--batches", "3",
                "--bootstrap", "12", "--output-dir", directory]; io=output)
            text = String(take!(output))

            @test status == 0
            @test isfile(joinpath(directory, "open_code_capacity_scan.csv"))
            @test isfile(joinpath(directory, "open_code_capacity_scan.pdf"))
            @test occursin("completed L=3", text)
            @test !isfile(joinpath(directory, "trajectory_scan.csv"))
        end
    end
end
