import Yao
using LinearAlgebra: dot
using Yao: apply!, nqubits, statevec, zero_state

function _operation_signature(operation)
    return (operation.gate, Tuple(operation.qubits))
end

_block_signature(block) = [
    _operation_signature(operation)
    for layer in block.layers for operation in layer.operations
]

_typ_site(code, x, y) = data_qubit_index(code, distance(code) - y + 1, x)

function _yao_operation_signature(block)
    if hasproperty(block, :ctrl_locs)
        return (:CNOT, (only(block.ctrl_locs), only(block.locs)))
    end
    gate = block.content === Yao.H ? :H : block.content === Yao.X ? :X : :Z
    return (gate, Tuple(block.locs))
end

function _pauli_expectation(state, support, pauli)
    mask = sum(UInt(1) << (qubit - 1) for qubit in support)
    if pauli === :Z
        return sum(abs2(amplitude) *
                   (isodd(count_ones(UInt(index - 1) & mask)) ? -1 : 1)
                   for (index, amplitude) in enumerate(state))
    end
    return sum(conj(state[index + 1]) * state[Int(xor(UInt(index), mask)) + 1]
               for index in 0:(length(state) - 1))
end

@testset "Rotated planar sequential encoder" begin
    required_api = (
        :PlaquetteEncoder,
        :rotated_planar_encoder,
        :plaquette_blocks,
        :gate_layers,
        :yao_encoder,
        :verify_encoder_tableau,
    )
    for name in required_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), required_api)
        @testset "constructor defaults, validation, and independent selectors" begin
            code = RotatedPlanarCode(3)
            default_encoder = rotated_planar_encoder(code)
            @test default_encoder.code === code
            @test default_encoder.logical_state === :zero
            @test default_encoder.construction === :bp
            for logical_state in (:zero, :one, :plus, :minus), construction in (:as, :bp)
                encoder = rotated_planar_encoder(
                    code; logical_state=logical_state, construction=construction)
                @test encoder.logical_state === logical_state
                @test encoder.construction === construction
            end
            @test_throws ArgumentError rotated_planar_encoder(code; logical_state=:cat)
            @test_throws ArgumentError rotated_planar_encoder(code; construction=:mixed)
        end

        @testset "family metadata, direction, and deterministic shell order" begin
            code = RotatedPlanarCode(3)
            bp = rotated_planar_encoder(code; construction=:bp)
            as = rotated_planar_encoder(code; construction=:as)
            bp_plaquettes = filter(block -> block.kind === :plaquette, plaquette_blocks(bp))
            as_plaquettes = filter(block -> block.kind === :plaquette, plaquette_blocks(as))

            @test !isempty(gate_layers(bp))
            @test !isempty(gate_layers(as))
            @test [_operation_signature(op) for layer in gate_layers(bp) for op in layer.operations] !=
                  [_operation_signature(op) for layer in gate_layers(as) for op in layer.operations]
            @test all(block.family === :bp && block.source_check === :B_p &&
                      block.direction === :outgoing_control for block in bp_plaquettes)
            @test all(block.family === :as && block.source_check === :A_s &&
                      block.direction === :incoming_control for block in as_plaquettes)
            @test count(block -> block.kind === :logical_sector,
                        plaquette_blocks(bp)) == 1
            @test count(block -> block.kind === :logical_sector,
                        plaquette_blocks(as)) >= 1

            # Literal execution order of the rotated d=3 template, expressed
            # in the default x_ns check enumeration.
            @test getproperty.(bp_plaquettes, :source_check_index) == [2, 3, 4, 1]
            @test getproperty.(as_plaquettes, :source_check_index) == [4, 1, 2, 3]
            @test only(filter(
                block -> block.source_check === :logical,
                plaquette_blocks(bp))).source_support == logical_x_support(code)
            @test only(filter(
                block -> block.source_check === :logical,
                plaquette_blocks(as))).source_support == logical_z_support(code)
            for d in (3, 5, 7), construction in (:as, :bp)
                code = RotatedPlanarCode(d)
                encoder = rotated_planar_encoder(code; construction=construction)
                repeat_encoder = rotated_planar_encoder(code; construction=construction)
                plaquettes = filter(
                    block -> block.kind === :plaquette, plaquette_blocks(encoder))
                @test getproperty.(plaquettes, :geometric_order) ==
                      collect(eachindex(plaquettes))
                @test getproperty.(plaquettes, :source_check_index) ==
                      getproperty.(filter(
                          block -> block.kind === :plaquette,
                          plaquette_blocks(repeat_encoder)), :source_check_index)
                for block in plaquettes
                    expected_checks = construction === :bp ?
                        b_p_checks(code) : a_s_checks(code)
                    @test block.source_support in expected_checks
                    @test block.representative in block.source_support
                    @test sort!(unique!([
                        qubit for layer in block.layers for operation in layer.operations
                        for qubit in operation.qubits])) == sort(block.source_support)
                end
                # In the B_p/X construction the H-control is a newly introduced
                # qubit.  The literal A_s/Z diagram instead reuses controls, so
                # imposing this B_p causal invariant on A_s would reject the
                # requested d=3 Typ circuit.
                if construction === :bp
                    touched = Set{Int}()
                    for block in plaquettes
                        @test !(block.representative in touched)
                        union!(touched, block.source_support)
                    end
                end
            end
        end

        @testset "all requested logical states retain source-check blocks" begin
            for d in (3, 5, 7), construction in (:as, :bp),
                logical_state in (:zero, :one, :plus, :minus)
                encoder = rotated_planar_encoder(
                    RotatedPlanarCode(d); construction=construction,
                    logical_state=logical_state)
                plaquettes = filter(
                    block -> block.kind === :plaquette, plaquette_blocks(encoder))
                expected_checks = construction === :bp ?
                    b_p_checks(encoder.code) : a_s_checks(encoder.code)
                @test all(block -> block.source_support in expected_checks &&
                                  block.representative in block.source_support,
                          plaquettes)
            end
        end

        @testset "d=3 Typ local growth gates" begin
            code = RotatedPlanarCode(3; boundary_orientation=:x_ew)
            q(x, y) = _typ_site(code, x, y)
            bp_blocks = filter(
                block -> block.kind === :plaquette,
                plaquette_blocks(rotated_planar_encoder(code; construction=:bp)))
            as_blocks = filter(
                block -> block.kind === :plaquette,
                plaquette_blocks(rotated_planar_encoder(code; construction=:as)))

            @test [_block_signature(block) for block in bp_blocks] == [
                [(:H, (q(2, 1),)), (:CNOT, (q(2, 1), q(1, 1))),
                 (:CNOT, (q(2, 1), q(2, 2))), (:CNOT, (q(2, 2), q(1, 2)))],
                [(:H, (q(1, 3),)), (:CNOT, (q(1, 3), q(1, 2)))],
                [(:H, (q(3, 1),)), (:CNOT, (q(3, 1), q(3, 2)))],
                [(:H, (q(3, 3),)), (:CNOT, (q(3, 3), q(3, 2))),
                 (:CNOT, (q(3, 3), q(2, 3))), (:CNOT, (q(2, 3), q(2, 2)))],
            ]
            @test [_block_signature(block) for block in as_blocks] == [
                [(:H, (q(2, 1),)), (:CNOT, (q(1, 1), q(2, 1)))],
                [(:H, (q(2, 2),)), (:CNOT, (q(1, 2), q(2, 2))),
                 (:CNOT, (q(1, 3), q(2, 2))), (:CNOT, (q(2, 3), q(2, 2)))],
                [(:H, (q(3, 1),)), (:CNOT, (q(2, 1), q(3, 1))),
                 (:CNOT, (q(2, 2), q(3, 1))), (:CNOT, (q(3, 2), q(3, 1)))],
                [(:H, (q(3, 3),)), (:CNOT, (q(2, 3), q(3, 3)))],
            ]
            @test getproperty.(bp_blocks, :representative) ==
                  [q(2, 1), q(1, 3), q(3, 1), q(3, 3)]
            @test getproperty.(as_blocks, :representative) ==
                  [q(2, 1), q(2, 2), q(3, 1), q(3, 3)]

            rotated = RotatedPlanarCode(3; boundary_orientation=:x_ns)
            rotate(qubit) = begin
                row, column = data_qubit_coordinate(code, qubit)
                data_qubit_index(rotated, column, distance(rotated) + 1 - row)
            end
            for construction in (:as, :bp)
                xew_blocks = filter(
                    block -> block.kind === :plaquette,
                    plaquette_blocks(rotated_planar_encoder(
                        code; construction=construction)))
                xns_blocks = filter(
                    block -> block.kind === :plaquette,
                    plaquette_blocks(rotated_planar_encoder(
                        rotated; construction=construction)))
                expected = [[(operation.gate, Tuple(rotate.(operation.qubits)))
                             for layer in block.layers for operation in layer.operations]
                            for block in xew_blocks]
                @test [_block_signature(block) for block in xns_blocks] == expected
            end
        end

        @testset "stored layers are exact, local, and bounded" begin
            for d in (3, 5, 7), construction in (:as, :bp)
                encoder = rotated_planar_encoder(
                    RotatedPlanarCode(d); construction=construction)
                @test [layer for block in plaquette_blocks(encoder) for layer in block.layers] ==
                      gate_layers(encoder)
                for layer in gate_layers(encoder)
                    @test layer.active_qubits ==
                          sort!(unique!(reduce(vcat, [copy(op.qubits) for op in layer.operations])))
                    for operation in layer.operations
                        @test all(qubit -> 1 <= qubit <= d^2, operation.qubits)
                        if operation.gate === :CNOT
                            @test length(operation.qubits) == 2
                            first_coordinate = data_qubit_coordinate(
                                encoder.code, operation.qubits[1])
                            second_coordinate = data_qubit_coordinate(
                                encoder.code, operation.qubits[2])
                            @test max(abs(first_coordinate[1] - second_coordinate[1]),
                                      abs(first_coordinate[2] - second_coordinate[2])) == 1
                        end
                    end
                end
            end
        end

        @testset "tableau verification covers distances and orientations" begin
            for d in (3, 5, 7), orientation in (:x_ns, :x_ew),
                construction in (:as, :bp), logical_state in (:zero, :one, :plus, :minus)
                code = RotatedPlanarCode(d; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=logical_state)
                @test verify_encoder_tableau(encoder)
            end
        end

        @testset "d=3 Yao states obey checks and logical sectors" begin
            for orientation in (:x_ns, :x_ew), logical_state in (:zero, :one, :plus, :minus)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                states = Dict{Symbol,Vector{ComplexF64}}()
                for construction in (:as, :bp)
                    encoder = rotated_planar_encoder(
                        code; construction=construction, logical_state=logical_state)
                    block = yao_encoder(encoder)
                    @test nqubits(block) == 9
                    stored_operations = [
                        _operation_signature(operation)
                        for layer in gate_layers(encoder) for operation in layer.operations
                    ]
                    @test _yao_operation_signature.(block.blocks) == stored_operations

                    register = zero_state(9)
                    apply!(register, block)
                    state = Vector{ComplexF64}(statevec(register))
                    states[construction] = state
                    for support in a_s_checks(code)
                        @test _pauli_expectation(state, support, :Z) ≈ 1 atol=1e-12
                    end
                    for support in b_p_checks(code)
                        @test _pauli_expectation(state, support, :X) ≈ 1 atol=1e-12
                    end
                    if logical_state in (:zero, :one)
                        expected = logical_state === :zero ? 1 : -1
                        @test _pauli_expectation(
                            state, logical_z_support(code), :Z) ≈ expected atol=1e-12
                    else
                        expected = logical_state === :plus ? 1 : -1
                        @test _pauli_expectation(
                            state, logical_x_support(code), :X) ≈ expected atol=1e-12
                    end
                end
                @test abs(dot(states[:as], states[:bp])) ≈ 1 atol=1e-12
            end
        end
    end
end
