using Test
using Random
using TopoNoise

toric3d_x_replay_test_syndrome(matrix, bits) =
    BitVector(isodd.(Int.(matrix) * Int.(bits)))

@testset "Paired diagonal X replay uses one immutable fault record" begin
    code = Toric3DEncodedSlice(2, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    slice_decoder = build_toric3d_slice_decoder(code.lattice)
    final_decoder = build_toric3d_x_threshold_decoder(code; p=0.02)
    matrix = plaquette_check_matrix(code)
    slice_matrix = toric3d_slice_check_matrix(code.lattice)
    clean = sample_toric3d_faults(MersenneTwister(7), encoder)

    for routing in (:one_step, :full_path)
        pair = toric3d_x_feedback_pair(encoder, clean, slice_decoder, final_decoder;
            routing, record_history=true)
        for result in (pair.final_only, pair.slice_feedback)
            @test !any(result.frame.x)
            @test !any(result.frame.z)
            @test !any(result.syndrome)
            @test !any(result.recovery)
            @test !any(result.residual)
            @test !result.logical_1
            @test !result.logical_2
            @test !result.failure
        end
    end

    for routing in (:one_step, :full_path),
        (gate, operation) in enumerate(encoder.operations), qubit in operation.qubits
        record = with_toric3d_fault(encoder, clean, gate, qubit, :X)
        saved_x, saved_z = copy(record.x), copy(record.z)
        pair = toric3d_x_feedback_pair(encoder, record, slice_decoder, final_decoder;
            routing, record_history=true)
        @test pair == toric3d_x_feedback_pair(encoder, record, slice_decoder,
            final_decoder; routing, record_history=true)
        @test record.x == saved_x
        @test record.z == saved_z
        for result in (pair.final_only, pair.slice_feedback)
            @test !any(result.frame.z)
            @test toric3d_x_replay_test_syndrome(matrix, result.recovery) ==
                result.syndrome
        end
        for checkpoint in pair.slice_feedback.history
            @test toric3d_x_replay_test_syndrome(slice_matrix,
                checkpoint.recovery) == checkpoint.syndrome
        end
    end

    odd = falses(size(slice_matrix, 1))
    odd[1] = true
    @test_throws ArgumentError decode_toric3d_slice(slice_decoder, odd)
    @test_throws ArgumentError replay_toric3d_x_feedback(encoder, clean,
        slice_decoder, final_decoder; feedback=true, routing=:invalid)

    other_slice = build_toric3d_slice_decoder(Toric3DSliceLattice(1, 2, 2))
    @test_throws DimensionMismatch replay_toric3d_x_feedback(encoder, clean,
        other_slice, final_decoder; feedback=true)
    other_final = build_toric3d_x_threshold_decoder(Toric3DEncodedSlice(1, 2, 2);
        p=0.02)
    @test_throws DimensionMismatch replay_toric3d_x_feedback(encoder, clean,
        slice_decoder, other_final; feedback=true)
    other_encoder = toric3d_slice_encoder(Toric3DSliceLattice(1, 2, 2))
    wrong_faults = sample_toric3d_faults(MersenneTwister(8), other_encoder)
    @test_throws ArgumentError replay_toric3d_x_feedback(encoder, wrong_faults,
        slice_decoder, final_decoder; feedback=true)
    z_fault = with_toric3d_fault(encoder, clean, 1,
        only(encoder.operations[1].qubits), :Z)
    @test_throws ArgumentError replay_toric3d_x_feedback(encoder, z_fault,
        slice_decoder, final_decoder; feedback=true)
end
