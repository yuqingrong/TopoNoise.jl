using Test
using Random
using TopoNoise

@testset "3D gate output faults" begin
    encoder = toric3d_encoder(ToricCode3D(2, 2, 2))
    n = data_qubit_count(encoder.code)
    clean = sample_toric3d_faults(MersenneTwister(1), encoder; p_x=0, p_z=0)
    @test size(clean.x) == size(clean.z) == (2, length(encoder.operations))
    @test !any(clean.x) && !any(clean.z)
    @test all(!clean.x[2, step] && !clean.z[2, step]
              for (step, op) in enumerate(encoder.operations) if op.gate === :H)
    @test sample_toric3d_faults(MersenneTwister(20), encoder;
        p_x=0.13, p_z=0.07) ==
          sample_toric3d_faults(MersenneTwister(20), encoder;
        p_x=0.13, p_z=0.07)
    @test_throws ArgumentError sample_toric3d_faults(MersenneTwister(1), encoder; p_x=-0.1)
    @test_throws ArgumentError with_toric3d_fault(encoder, clean, 1, 1, :X)
    @test_throws ArgumentError with_toric3d_fault(encoder, clean, 1,
        encoder.operations[1].qubits[1], :Q)
    @test_throws ArgumentError Toric3DFaultRecord(0, 2, 2, clean.x, clean.z)
    @test_throws ArgumentError Toric3DFaultRecord(2, 1, 2, clean.x, clean.z)
    @test_throws ArgumentError Toric3DFaultRecord(2, 2, 1, clean.x, clean.z)
    @test_throws ArgumentError Toric3DFaultRecord(2, 2, 2,
        falses(1, length(encoder.operations)), clean.z)
    for (step, operation) in enumerate(encoder.operations)
        for qubit in operation.qubits, pauli in (:X, :Z)
            record = with_toric3d_fault(encoder, clean, step, qubit, pauli)
            frame = PauliFrame(falses(n), falses(n))
            propagate_toric3d_segment!(frame, encoder, record, 1, step)
            @test count(frame.x) + count(frame.z) == 1
            @test (pauli === :X ? frame.x : frame.z)[qubit]
            propagate_toric3d_segment!(frame, encoder, record,
                step+1, length(encoder.operations))
            # An independently evolved basis vector checks every eligible fault.
            x, z = falses(n), falses(n)
            (pauli === :X ? x : z)[qubit] = true
            for later in encoder.operations[step+1:end]
                if later.gate === :H
                    q = only(later.qubits)
                    x[q], z[q] = z[q], x[q]
                else
                    c, t = later.qubits
                    x[t] ⊻= x[c]
                    z[c] ⊻= z[t]
                end
            end
            @test frame == PauliFrame(x, z)
            syndrome = toric3d_syndrome(encoder.code, frame)
            @test syndrome.stars == BitVector(isodd(count(frame.z[s])) for s in star_checks(encoder.code))
            @test syndrome.plaquettes == BitVector(isodd(count(frame.x[s])) for s in plaquette_checks(encoder.code))
            parities = toric3d_logical_parities(encoder.code, frame)
            @test parities.logical_x == isodd(count(frame.x[logical_z_support(encoder.code)]))
            @test parities.logical_z == isodd(count(frame.z[logical_x_support(encoder.code)]))
        end
    end
    first_layer_end = encoder.checkpoints[1]
    crossing_qubit = encoder.operations[first_layer_end].qubits[1]
    crossing = with_toric3d_fault(encoder, clean, first_layer_end,
                                  crossing_qubit, :Z)
    frame = PauliFrame(falses(n), falses(n))
    propagate_toric3d_segment!(frame, encoder, crossing, 1, first_layer_end)
    @test frame.z[crossing_qubit]
    later_control = data_qubit_index(encoder.code, :x, 2, 1, 1)
    @test !frame.z[later_control]
    propagate_toric3d_segment!(frame, encoder, crossing,
        first_layer_end + 1, length(encoder.operations))
    @test frame.z[later_control]
    @test !any(clean.x) && !any(clean.z)
    y_record = with_toric3d_fault(encoder, clean, 1,
        encoder.operations[1].qubits[1], :Y)
    @test y_record.x[1, 1] && y_record.z[1, 1]
    @test with_toric3d_fault(encoder, y_record, 1,
        encoder.operations[1].qubits[1], :Y) == clean
    bad_h = Toric3DFaultRecord(2, 2, 2, clean.x, clean.z)
    bad_h.x[2, 1] = true
    @test_throws ArgumentError propagate_toric3d_segment!(
        PauliFrame(falses(n), falses(n)), encoder, bad_h, 1, 1)
end
