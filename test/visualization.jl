using CairoMakie: CairoMakie, Figure
import Yao

function rendered_text(fig)
    labels = String[]
    for block in fig.content
        hasproperty(block, :scene) || continue
        for plot in block.scene.plots
            hasproperty(plot, :text) || continue
            append!(labels, string.(plot.text[]))
        end
    end
    return labels
end

@testset "Native Yao circuit figures" begin
    @test isdefined(TopoNoise, :yao_circuit)

    if isdefined(TopoNoise, :yao_circuit)
        mktempdir() do directory
            for (rows, cols) in ((1, 1), (2, 2))
                block = yao_circuit(toric_code_peps(rows, cols))
                @test Yao.plot(block) !== nothing

                for format in (:svg, :pdf)
                    path = joinpath(directory, "yao-$(rows)x$(cols).$(format)")
                    Yao.vizcircuit(block; format=format, filename=path)
                    @test isfile(path)
                    @test filesize(path) > 100
                end
            end
            @test startswith(
                read(joinpath(directory, "yao-1x1.svg"), String), "<?xml")
            @test startswith(
                read(joinpath(directory, "yao-2x2.pdf"), String), "%PDF")
        end
    end
end

@testset "Toric-code vector figures" begin
    @test isdefined(TopoNoise, :plot_peps_graph)
    @test isdefined(TopoNoise, :plot_sequential_circuit)

    if isdefined(TopoNoise, :plot_peps_graph)
        peps = toric_code_peps(2, 2)
        circuit = sequential_circuit_graph(peps)
        peps_figure = plot_peps_graph(peps; show_index_labels=true)
        circuit_figure = plot_sequential_circuit(circuit)
        expanded_figure = plot_sequential_circuit(
            circuit; show_layer_labels=false, expand_physical_buses=true)

        @test peps_figure isa Figure
        @test circuit_figure isa Figure
        @test expanded_figure isa Figure

        peps_labels = rendered_text(peps_figure)
        @test all("T[$row,$col]" in peps_labels for row in 1:2 for col in 1:2)
        @test all(any(startswith("$direction:"), peps_labels)
                  for direction in ("E", "N", "W", "S"))

        circuit_labels = rendered_text(circuit_figure)
        @test all("U[$row,$col]" in circuit_labels for row in 1:2 for col in 1:2)
        @test any(contains("W→E"), circuit_labels)
        @test any(contains("S→N"), circuit_labels)
        @test "|0⟩" in circuit_labels
        @test "|+⟩" in circuit_labels
        @test "⟨+|" in circuit_labels
        @test all("layer $layer" in circuit_labels for layer in 1:3)

        mktempdir() do directory
            outputs = (
                joinpath(directory, "peps.svg") => peps_figure,
                joinpath(directory, "peps.pdf") => peps_figure,
                joinpath(directory, "circuit.svg") => circuit_figure,
                joinpath(directory, "circuit.pdf") => circuit_figure,
            )
            for (path, figure) in outputs
                CairoMakie.save(path, figure)
                @test isfile(path)
                @test filesize(path) > 100
            end
            @test startswith(read(joinpath(directory, "peps.svg"), String), "<?xml")
            @test startswith(read(joinpath(directory, "circuit.pdf"), String), "%PDF")
        end
    end
end
