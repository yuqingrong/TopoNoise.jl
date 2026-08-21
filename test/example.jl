const _EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "examples", "generate_toric_circuit.jl")

@testset "Runnable toric-code circuit example" begin
    @test isfile(_EXAMPLE_PATH)

    if isfile(_EXAMPLE_PATH)
        include(_EXAMPLE_PATH)

        default_output = GenerateToricCircuitExample.generate(1, 1).outputs
        @test all(startswith(path, normpath(joinpath(
            @__DIR__, "..", "results", "toric-circuit"))) for path in default_output)

        mktempdir() do directory
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                String[]; output_dir=directory, io=output, error_io=output) == 0

            text = String(take!(output))
            @test contains(text, "lattice: 2x2")
            @test contains(text, "qubits: 20")
            @test contains(text, "layer 1: (2,1)")
            @test contains(text, "layer 2: (2,2), (1,1)")
            @test contains(text, "layer 3: (1,2)")

            for extension in ("svg", "pdf", "png")
                path = joinpath(directory, "toric_circuit_2x2.$extension")
                @test isfile(path)
                @test filesize(path) > 100
                @test contains(text, "wrote: $(abspath(path))")
            end
            @test startswith(
                read(joinpath(directory, "toric_circuit_2x2.svg"), String), "<?xml")
            @test startswith(
                read(joinpath(directory, "toric_circuit_2x2.pdf"), String), "%PDF")
        end

        mktempdir() do directory
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                ["1", "1"];
                output_dir=directory, io=output, error_io=output) == 0
            @test contains(String(take!(output)), "qubits: 6")
            @test isfile(joinpath(directory, "toric_circuit_1x1.svg"))
        end

        for arguments in (
                ["2"], ["1", "2", "3"], ["two", "2"],
                ["0", "2"], ["2", "-1"])
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                arguments;
                output_dir="unused", io=output, error_io=output) == 1
            text = String(take!(output))
            @test contains(text, "Usage:")
            @test contains(text, "positive integers")
        end
    end
end
