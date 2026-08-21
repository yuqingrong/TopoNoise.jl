using Test

const _ROTATED_CODE_CAPACITY_EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "examples", "scan_rotated_code_capacity.jl")

@testset "rotated code-capacity example" begin
    @test isfile(_ROTATED_CODE_CAPACITY_EXAMPLE_PATH)

    if isfile(_ROTATED_CODE_CAPACITY_EXAMPLE_PATH)
        include(_ROTATED_CODE_CAPACITY_EXAMPLE_PATH)

        defaults = ScanRotatedCodeCapacityExample._parse_options([
            "--p-min", "0.05", "--p-max", "0.10", "--p-step", "0.05"])
        @test defaults.output_dir == normpath(joinpath(
            @__DIR__, "..", "results", "rotated-code-capacity"))

        mktempdir() do directory
            output = IOBuffer()
            status = ScanRotatedCodeCapacityExample.main([
                "--p-min", "0.05", "--p-max", "0.10", "--p-step", "0.05",
                "--sizes", "3,5", "--shots", "12", "--batches", "3",
                "--bootstrap", "12", "--seed", "31", "--output-dir", directory];
                io=output)

            @test status == 0
            @test all(isfile(joinpath(directory, name)) for name in (
                "rotated_code_capacity_scan.csv",
                "rotated_code_capacity_crossings.csv",
                "rotated_code_capacity_scan.svg",
                "rotated_code_capacity_scan.pdf",
                "rotated_code_capacity_scan.png"))
            @test occursin("unbracketed", String(take!(output)))
            @test startswith(read(joinpath(
                directory, "rotated_code_capacity_scan.svg"), String), "<?xml")
            @test startswith(read(joinpath(
                directory, "rotated_code_capacity_scan.pdf"), String), "%PDF")
            @test read(joinpath(directory, "rotated_code_capacity_scan.png"))[1:8] ==
                UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
        end
    end
end
