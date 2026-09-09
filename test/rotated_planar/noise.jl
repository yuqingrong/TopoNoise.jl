@testset "Rotated planar circuit Pauli noise" begin
    noise_api = (
        :CircuitPauliNoise,
        :CircuitFaultStep,
        :CircuitFaultRecord,
        :sample_fault_record,
    )
    for name in noise_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), noise_api)
        @testset "noise rates, defaults, and validation" begin
            default = CircuitPauliNoise(0.125)
            @test default.p_x == 0.125
            @test default.p_z == 0.125
            @test default.clock === :gate_layer

            unequal = CircuitPauliNoise(0.1; p_x=0.2, p_z=0.3, clock=:plaquette)
            @test unequal.p_x == 0.2
            @test unequal.p_z == 0.3
            @test unequal.clock === :plaquette

            for invalid in (-0.1, 0.5, Inf, -Inf, NaN)
                @test_throws ArgumentError CircuitPauliNoise(0.1; p_x=invalid)
                @test_throws ArgumentError CircuitPauliNoise(0.1; p_z=invalid)
            end
            @test_throws ArgumentError CircuitPauliNoise(0.1; clock=:round)
            @test_throws ArgumentError CircuitPauliNoise(-0.1, 0.2, :gate_layer)
            @test_throws ArgumentError CircuitPauliNoise(0.1, 0.2, :round)
        end

        @testset "records retain exact eligible-step identity and support" begin
            for construction in (:as, :bp), clock in (:gate_layer, :plaquette)
                code = RotatedPlanarCode(3)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=:minus)
                record = sample_fault_record(
                    MersenneTwister(11), encoder, CircuitPauliNoise(0; clock=clock))

                @test record.clock === clock
                @test record.distance == 3
                @test record.boundary_orientation === :x_ns
                @test record.construction === construction
                @test record.logical_state === :minus
                expected = clock === :gate_layer ?
                    [(i, layer.block_kind, layer.block_order, layer.active_qubits)
                     for (i, layer) in enumerate(gate_layers(encoder))] :
                    [(i, block.kind, block.geometric_order, block.source_support)
                     for (i, block) in enumerate(filter(
                         candidate -> candidate.kind === :plaquette,
                         plaquette_blocks(encoder)))]
                actual = [(step.step_index, step.block_kind, step.block_order,
                           collect(step.support))
                          for step in record.steps]
                @test actual == expected
                @test all(step -> length(step.x) == length(step.support) &&
                                  length(step.z) == length(step.support), record.steps)
                @test all(step -> !any(step.x) && !any(step.z), record.steps)
            end
        end

        @testset "sampling is reproducible and X/Z paths are independent" begin
            encoder = rotated_planar_encoder(RotatedPlanarCode(3); construction=:bp)
            for clock in (:gate_layer, :plaquette)
                noise = CircuitPauliNoise(0.1; p_x=0.17, p_z=0.31, clock=clock)
                @test sample_fault_record(MersenneTwister(90210), encoder, noise) ==
                      sample_fault_record(MersenneTwister(90210), encoder, noise)

                only_z = sample_fault_record(
                    MersenneTwister(7), encoder,
                    CircuitPauliNoise(0; p_x=0, p_z=0.49, clock=clock))
                @test !any(step -> any(step.x), only_z.steps)
                @test any(step -> any(step.z), only_z.steps)

                only_x = sample_fault_record(
                    MersenneTwister(7), encoder,
                    CircuitPauliNoise(0; p_x=0.49, p_z=0, clock=clock))
                @test any(step -> any(step.x), only_x.steps)
                @test !any(step -> any(step.z), only_x.steps)
            end
        end

        @testset "coincident draws remain a Y event" begin
            step = CircuitFaultStep(1, :plaquette, 1, [5], trues(1), trues(1))
            @test step.x[1]
            @test step.z[1]
        end
    end
end

function _literal_encoder(operations)
    code = RotatedPlanarCode(3)
    layers = [
        TopoNoise.EncoderGateLayer(
            :logical_sector, 0, [operation], sort(unique(operation.qubits)))
        for operation in operations
    ]
    block = TopoNoise.PlaquetteEncoderBlock(
        :logical_sector, :bp, :logical, 0, Int[], Int[], Int[], 1,
        (0, 0), 0, (typemin(Int), 0, 0, 0), 0, :outgoing_control, layers)
    return PlaquetteEncoder(code, :zero, :bp, [block], layers)
end

function _record_with_fault(encoder, step_index, pauli)
    record = sample_fault_record(
        MersenneTwister(1), encoder, CircuitPauliNoise(0; clock=:gate_layer))
    qubit = record.steps[step_index].support[1]
    return with_pauli_fault(record, step_index, qubit, pauli)
end

function _with_test_fault(record, step_index, pauli)
    qubit = record.steps[step_index].support[1]
    return with_pauli_fault(record, step_index, qubit, pauli)
end

@testset "Fault records are authoritative immutable values" begin
    @test isdefined(TopoNoise, :with_pauli_fault)

    encoder = rotated_planar_encoder(RotatedPlanarCode(3))
    sampled = sample_fault_record(
        MersenneTwister(70), encoder, CircuitPauliNoise(0))
    source_step = sampled.steps[end]

    support_input = collect(source_step.support)
    x_input = collect(source_step.x)
    z_input = collect(source_step.z)
    original_support = Tuple(support_input)
    step = CircuitFaultStep(
        source_step.step_index, source_step.block_kind, source_step.block_order,
        support_input, x_input, z_input)
    support_input[1] = support_input[1] == 1 ? 2 : 1
    x_input[1] = true
    z_input[1] = true
    @test Tuple(step.support) == original_support
    @test !any(step.x)
    @test !any(step.z)

    steps_input = collect(sampled.steps)
    record = CircuitFaultRecord(
        sampled.clock, sampled.distance, sampled.boundary_orientation,
        sampled.construction, sampled.logical_state, steps_input)
    replay_before = propagate_pauli_frame(encoder, record)
    steps_input[end] = CircuitFaultStep(
        source_step.step_index, source_step.block_kind, source_step.block_order,
        collect(source_step.support), trues(length(source_step.support)),
        falses(length(source_step.support)))
    @test propagate_pauli_frame(encoder, record) == replay_before

    @test record.steps isa Tuple
    @test record.steps[end].support isa Tuple
    @test record.steps[end].x isa Tuple
    @test record.steps[end].z isa Tuple
    @test_throws MethodError setindex!(record.steps, source_step, 1)
    @test_throws MethodError setindex!(record.steps[end].x, true, 1)

    qubit = record.steps[end].support[1]
    x_record = with_pauli_fault(record, length(record.steps), qubit, :X)
    y_record = with_pauli_fault(record, length(record.steps), qubit, :Y)
    @test !record.steps[end].x[1] && !record.steps[end].z[1]
    @test x_record.steps[end].x[1] && !x_record.steps[end].z[1]
    @test y_record.steps[end].x[1] && y_record.steps[end].z[1]
    @test_throws ArgumentError with_pauli_fault(
        record, length(record.steps), qubit, :W)
    @test_throws ArgumentError with_pauli_fault(record, 0, qubit, :X)
    @test_throws ArgumentError with_pauli_fault(
        record, length(record.steps), 99, :X)
end

@testset "Rotated planar Pauli propagation and syndrome" begin
    frame_api = (
        :PauliFrame,
        :SyndromeRecord,
        :propagate_pauli_frame,
        :measure_syndrome,
        :syndrome_bits,
    )
    for name in frame_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), frame_api)
        @testset "literal H swaps frame axes and ideal X does not" begin
            encoder = _literal_encoder([
                TopoNoise.EncoderOperation(:X, [1]),
                TopoNoise.EncoderOperation(:H, [1]),
            ])
            @test propagate_pauli_frame(
                encoder, _record_with_fault(encoder, 1, :X)) ==
                  PauliFrame(falses(9), BitVector([true; falses(8)]))
            @test propagate_pauli_frame(
                encoder, _record_with_fault(encoder, 1, :Z)) ==
                  PauliFrame(BitVector([true; falses(8)]), falses(9))
            @test propagate_pauli_frame(
                encoder, _record_with_fault(encoder, 1, :Y)) ==
                  PauliFrame(BitVector([true; falses(8)]),
                             BitVector([true; falses(8)]))
        end

        @testset "literal CNOT propagates each binary direction" begin
            control_fixture = _literal_encoder([
                TopoNoise.EncoderOperation(:X, [1]),
                TopoNoise.EncoderOperation(:CNOT, [1, 2]),
            ])
            target_fixture = _literal_encoder([
                TopoNoise.EncoderOperation(:X, [2]),
                TopoNoise.EncoderOperation(:CNOT, [1, 2]),
            ])
            x_control = propagate_pauli_frame(
                control_fixture, _record_with_fault(control_fixture, 1, :X))
            z_control = propagate_pauli_frame(
                control_fixture, _record_with_fault(control_fixture, 1, :Z))
            x_target = propagate_pauli_frame(
                target_fixture, _record_with_fault(target_fixture, 1, :X))
            z_target = propagate_pauli_frame(
                target_fixture, _record_with_fault(target_fixture, 1, :Z))
            @test findall(x_control.x) == [1, 2]
            @test findall(z_control.z) == [1]
            @test findall(x_target.x) == [2]
            @test findall(z_target.z) == [1, 2]
        end

        @testset "output X, Z, and Y address the CSS-dual checks" begin
            code = RotatedPlanarCode(3)
            x = falses(9)
            z = falses(9)
            x[5] = true
            x_syndrome = measure_syndrome(code, PauliFrame(x, z))
            @test x_syndrome.a_s == BitVector([true, true, false, false])
            @test x_syndrome.b_p == falses(4)

            x[5] = false
            z[5] = true
            z_syndrome = measure_syndrome(code, PauliFrame(x, z))
            @test z_syndrome.a_s == falses(4)
            @test z_syndrome.b_p == BitVector([true, true, false, false])

            x[5] = true
            y_syndrome = measure_syndrome(code, PauliFrame(x, z))
            @test y_syndrome.a_s == BitVector([true, true, false, false])
            @test y_syndrome.b_p == BitVector([true, true, false, false])
            @test syndrome_bits(y_syndrome) ==
                  BitVector([true, true, false, false,
                             true, true, false, false])
        end

        @testset "noiseless schedules and coincident output faults" begin
            for orientation in (:x_ns, :x_ew), construction in (:as, :bp),
                logical_state in (:zero, :one, :plus, :minus),
                clock in (:gate_layer, :plaquette)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=logical_state)
                record = sample_fault_record(
                    MersenneTwister(2), encoder, CircuitPauliNoise(0; clock=clock))
                frame = propagate_pauli_frame(encoder, record)
                @test !any(frame.x)
                @test !any(frame.z)
                @test !any(syndrome_bits(measure_syndrome(code, frame)))
            end

            encoder = rotated_planar_encoder(RotatedPlanarCode(3))
            record = sample_fault_record(
                MersenneTwister(3), encoder,
                CircuitPauliNoise(0; clock=:gate_layer))
            final_step = record.steps[end]
            x = falses(length(final_step.support))
            z = falses(length(final_step.support))
            x[1] = true
            z[1] = true
            coincident_step = CircuitFaultStep(
                final_step.step_index, final_step.block_kind,
                final_step.block_order, collect(final_step.support), x, z)
            manual_steps = collect(record.steps)
            manual_steps[end] = coincident_step
            manual_record = CircuitFaultRecord(
                record.clock, record.distance, record.boundary_orientation,
                record.construction, record.logical_state, manual_steps)
            qubit = coincident_step.support[1]
            frame = propagate_pauli_frame(encoder, manual_record)
            @test frame.x[qubit] && frame.z[qubit]
            @test count(frame.x) == 1
            @test count(frame.z) == 1
        end

        @testset "record compatibility is validated before replay" begin
            encoder = rotated_planar_encoder(RotatedPlanarCode(3))
            record = sample_fault_record(
                MersenneTwister(4), encoder, CircuitPauliNoise(0))
            wrong_distance = CircuitFaultRecord(
                record.clock, 5, record.boundary_orientation,
                record.construction, record.logical_state, record.steps)
            @test_throws ArgumentError propagate_pauli_frame(encoder, wrong_distance)

            @test_throws ArgumentError CircuitFaultStep(
                record.steps[1].step_index, record.steps[1].block_kind,
                record.steps[1].block_order, record.steps[1].support,
                falses(length(record.steps[1].support) + 1), record.steps[1].z)
            @test_throws ArgumentError measure_syndrome(
                encoder.code, PauliFrame(falses(8), falses(8)))
        end
    end
end

@testset "Rotated planar Yao ancilla syndrome oracle" begin
    @test isdefined(TopoNoise, :sample_yao_syndrome)

    if isdefined(TopoNoise, :sample_yao_syndrome)
        @testset "manual output X, Z, and Y faults match algebra" begin
            for orientation in (:x_ns, :x_ew), construction in (:as, :bp),
                logical_state in (:zero, :plus), clock in (:gate_layer, :plaquette),
                pauli in (:X, :Z, :Y)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=logical_state)
                record = sample_fault_record(
                    MersenneTwister(20), encoder, CircuitPauliNoise(0; clock=clock))
                record = _with_test_fault(record, length(record.steps), pauli)
                algebraic = measure_syndrome(
                    code, propagate_pauli_frame(encoder, record))
                @test sample_yao_syndrome(
                    MersenneTwister(30), code, encoder, record) == algebraic
            end
        end

        @testset "seeded random records match algebra without resampling" begin
            for orientation in (:x_ns, :x_ew), construction in (:as, :bp),
                logical_state in (:zero, :plus), clock in (:gate_layer, :plaquette),
                seed in (101, 202)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=logical_state)
                record = sample_fault_record(
                    MersenneTwister(seed), encoder,
                    CircuitPauliNoise(0; p_x=0.13, p_z=0.29, clock=clock))
                algebraic = measure_syndrome(
                    code, propagate_pauli_frame(encoder, record))
                @test sample_yao_syndrome(
                    MersenneTwister(seed + 1), code, encoder, record) == algebraic
                @test sample_yao_syndrome(
                    MersenneTwister(seed + 2), code, encoder, record) == algebraic
            end
        end

        @testset "noiseless extraction is zero and larger codes are rejected" begin
            for orientation in (:x_ns, :x_ew), construction in (:as, :bp),
                logical_state in (:zero, :one, :plus, :minus),
                clock in (:gate_layer, :plaquette)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(
                    code; construction=construction, logical_state=logical_state)
                record = sample_fault_record(
                    MersenneTwister(40), encoder, CircuitPauliNoise(0; clock=clock))
                @test !any(syndrome_bits(sample_yao_syndrome(
                    MersenneTwister(41), code, encoder, record)))
            end

            code = RotatedPlanarCode(5)
            encoder = rotated_planar_encoder(code)
            record = sample_fault_record(
                MersenneTwister(50), encoder, CircuitPauliNoise(0))
            @test_throws ArgumentError sample_yao_syndrome(
                MersenneTwister(51), code, encoder, record)
        end
    end
end
