import Yao
using Yao: H, apply!, nqubits, product_state, repeat, statevec, zero_state

@testset "Toric-code sequential circuit schedule" begin
    @test isdefined(TopoNoise, :ToricCodeSequentialCircuit)
    @test isdefined(TopoNoise, :sequential_circuit_graph)
    @test isdefined(TopoNoise, :circuit_layers)
    @test isdefined(TopoNoise, :yao_unitary)
    @test isdefined(TopoNoise, :yao_circuit)
    @test isdefined(TopoNoise, :log2_postselection_probability)
    @test isdefined(TopoNoise, :postselection_probability)

    if isdefined(TopoNoise, :ToricCodeSequentialCircuit)
        for (rows, cols) in ((1, 1), (1, 2), (2, 1), (2, 2), (3, 3))
            peps = toric_code_peps(rows, cols)
            circuit = sequential_circuit_graph(peps)
            layers = circuit_layers(circuit)
            vertices = rows * cols
            edges = rows * (cols - 1) + (rows - 1) * cols

            @test length(layers) == rows + cols - 1
            @test Set(vcat(layers...)) ==
                  Set((row, col) for row in 1:rows for col in 1:cols)
            @test all(length(unique(first.(layer))) == length(layer) for layer in layers)
            @test all(length(unique(last.(layer))) == length(layer) for layer in layers)
            @test circuit.horizontal_wires ==
                  (4 * vertices + 1):(4 * vertices + rows)
            @test circuit.vertical_wires ==
                  (4 * vertices + rows + 1):(4 * vertices + rows + cols)
            gates = Dict(gate.site => gate for layer in circuit.layers for gate in layer)
            for row in 1:rows, col in 1:cols
                first_physical = 4 * ((row - 1) * cols + col - 1) + 1
                @test gates[(row, col)].wires == (
                    first_physical, first_physical + 1,
                    first_physical + 2, first_physical + 3,
                    4 * vertices + row, 4 * vertices + rows + col)
            end
            scheduled_gates = [gate for layer in circuit.layers for gate in layer]
            embedded_blocks = yao_unitary(circuit).blocks
            @test length(embedded_blocks) == vertices
            for (block, gate) in zip(embedded_blocks, scheduled_gates)
                @test block.locs == gate.wires
                @test block.content.tag == "U[$(gate.site[1]),$(gate.site[2])]"
            end
            @test nqubits(yao_unitary(circuit)) == 4 * vertices + rows + cols
            @test circuit.carrier_input_state == :plus
            @test circuit.carrier_output_effect == :plus_projection
            @test log2_postselection_probability(circuit) == edges - 2 * vertices
            @test postselection_probability(circuit) == exp2(edges - 2 * vertices)
        end

        circuit = sequential_circuit_graph(toric_code_peps(3, 3))
        @test circuit_layers(circuit) == [
            [(3, 1)],
            [(3, 2), (2, 1)],
            [(3, 3), (2, 2), (1, 1)],
            [(2, 3), (1, 2)],
            [(1, 3)],
        ]
    end
end

@testset "Executable Yao circuit" begin
    @test isdefined(TopoNoise, :yao_circuit)

    if isdefined(TopoNoise, :yao_circuit)
        for (rows, cols) in ((1, 1), (1, 2), (2, 1))
            peps = toric_code_peps(rows, cols; scalar_type=ComplexF64)
            circuit = sequential_circuit_graph(peps)
            from_peps = yao_circuit(peps)
            from_circuit = yao_circuit(circuit)
            vertices = rows * cols
            wire_count = 4 * vertices + rows + cols
            carriers = (4 * vertices + 1):wire_count

            @test from_peps isa Yao.ChainBlock
            @test from_circuit isa Yao.ChainBlock
            @test nqubits(from_peps) == wire_count
            @test length(from_circuit.blocks) == 2
            @test from_circuit.blocks[1].content == H
            @test from_circuit.blocks[1].locs == Tuple(carriers)

            direct_register = zero_state(wire_count)
            apply!(direct_register, from_peps)

            manual_register = zero_state(wire_count)
            apply!(manual_register, repeat(wire_count, H, carriers))
            apply!(manual_register, yao_unitary(circuit))
            @test statevec(direct_register) ≈ statevec(manual_register)

            model_register = zero_state(wire_count)
            apply!(model_register, from_circuit)
            @test statevec(direct_register) ≈ statevec(model_register)
        end

        block = yao_circuit(toric_code_peps(3, 3))
        core = block.blocks[2]
        @test [gate.content.tag for gate in core.blocks] == [
            "U[3,1]",
            "U[3,2]", "U[2,1]",
            "U[3,3]", "U[2,2]", "U[1,1]",
            "U[2,3]", "U[1,2]",
            "U[1,3]",
        ]
    end
end

@testset "Yao local-gate convention and small-patch states" begin
    if isdefined(TopoNoise, :ToricCodeSequentialCircuit)
        local_matrix = reshape(toric_code_local_gate(ComplexF64), 64, 64)
        local_block = TopoNoise._yao_local_gate(ComplexF64, (1, 1))
        for input in (0, 1, 3, 17, 63)
            register = product_state(6, input)
            apply!(register, local_block)
            @test statevec(register) ≈ local_matrix[:, input + 1]
        end

        for (rows, cols) in ((1, 1), (1, 2), (2, 1))
            peps = toric_code_peps(rows, cols; scalar_type=ComplexF64)
            circuit = sequential_circuit_graph(peps)
            vertices = rows * cols
            wire_count = 4 * vertices + rows + cols
            carriers = (4 * vertices + 1):wire_count

            register = zero_state(wire_count)
            apply!(register, repeat(wire_count, H, carriers))
            apply!(register, yao_unitary(circuit))
            apply!(register, repeat(wire_count, H, carriers))

            projected = reshape(statevec(register), ntuple(_ -> 2, wire_count))[
                ntuple(wire -> wire <= 4 * vertices ? Colon() : 1, wire_count)...]
            circuit_state = vec(projected)
            @test sum(abs2, circuit_state) ≈ postselection_probability(circuit)
            circuit_state ./= sqrt(sum(abs2, circuit_state))

            physical = reduce(vcat, collect.(vec(peps.physical_indices)))
            peps_tensor = ITensors.contract(collect(vec(peps.tensors)))
            peps_state = vec(Array(peps_tensor, physical...))
            peps_state ./= sqrt(sum(abs2, peps_state))
            @test circuit_state ≈ peps_state atol=1e-11 rtol=1e-11
        end
    end
end
