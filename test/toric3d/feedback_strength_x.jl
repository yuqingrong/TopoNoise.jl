using Test
using Random
using TopoNoise

toric3d_strength_x_syndrome(matrix, bits) =
    BitVector(isodd.(Int.(matrix) * Int.(bits)))

@testset "X strength baselines preserve the established replay exactly" begin
    code = Toric3DEncodedSlice(2, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_x_feedback_decoders(code; p=0.02)
    @test decoders isa Toric3DXFeedbackDecoderFamily
    @test decoders.code === code
    @test decoders.slice.code === code.lattice
    @test length(decoders.prefixes) == 2
    @test decoders.prefixes[end].code.lattice.steps == 2

    clean = sample_toric3d_faults(MersenneTwister(7), encoder)
    random_record = sample_toric3d_faults(MersenneTwister(19), encoder;
        p_x=0.17, p_z=0.0)
    records = Toric3DFaultRecord[clean, random_record]
    append!(records,
        [with_toric3d_fault(encoder, clean, gate, qubit, :X)
         for (gate, operation) in enumerate(encoder.operations)
         for qubit in operation.qubits])

    for routing in (:one_step, :full_path), record in records
        established = toric3d_x_feedback_pair(encoder, record, decoders.slice,
            decoders.prefixes[end]; routing, record_history=true)
        @test replay_toric3d_x_strength(encoder, record, decoders;
            protocol=:none, routing, record_history=true) == established.final_only
        @test replay_toric3d_x_strength(encoder, record, decoders;
            protocol=:slice_local, routing, record_history=true) == established.slice_feedback
    end

    saved_x, saved_z = copy(random_record.x), copy(random_record.z)
    triplet = toric3d_x_strength_triplet(encoder, random_record, decoders;
        routing=:one_step, record_history=true)
    @test keys(triplet) == (:none, :slice_local, :full_prefix)
    @test random_record.x == saved_x
    @test random_record.z == saved_z
end

@testset "Full-prefix X feedback validates every checkpoint and final recovery" begin
    code = Toric3DEncodedSlice(2, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_x_feedback_decoders(code; p=0.02)
    clean = sample_toric3d_faults(MersenneTwister(23), encoder)

    clean_result = replay_toric3d_x_strength(encoder, clean, decoders;
        protocol=:full_prefix, record_history=true)
    @test !any(clean_result.frame.x)
    @test !any(clean_result.frame.z)
    @test !any(clean_result.syndrome)
    @test !any(clean_result.recovery)
    @test !any(clean_result.residual)
    @test !clean_result.failure
    @test length(clean_result.history) == 1

    random_record = sample_toric3d_faults(MersenneTwister(29), encoder;
        p_x=0.2, p_z=0.0)
    random_result = replay_toric3d_x_strength(encoder, random_record, decoders;
        protocol=:full_prefix, record_history=true)
    @test length(random_result.history) == 1
    for checkpoint in random_result.history
        decoder = decoders.prefixes[checkpoint.checkpoint]
        @test toric3d_strength_x_syndrome(decoder.matrix, checkpoint.recovery) ==
              checkpoint.pre_syndrome
        @test length(checkpoint.qubits) == size(decoder.matrix, 2)
        @test length(unique(checkpoint.qubits)) == length(checkpoint.qubits)
        @test !any(checkpoint.post_syndrome)
        @test !any(checkpoint.frontier_syndrome)
        @test checkpoint.future_before == checkpoint.future_after
        @test all(data_qubit_coordinate(code, q)[2] > checkpoint.checkpoint
                  for q in checkpoint.future_qubits)
    end

    for (gate, operation) in enumerate(encoder.operations), qubit in operation.qubits
        record = with_toric3d_fault(encoder, clean, gate, qubit, :X)
        result = replay_toric3d_x_strength(encoder, record, decoders;
            protocol=:full_prefix)
        final_decoder = decoders.prefixes[end]
        @test toric3d_strength_x_syndrome(final_decoder.matrix, result.recovery) ==
              result.syndrome
        @test !any(toric3d_strength_x_syndrome(final_decoder.matrix, result.residual))
        @test !any(result.frame.z)
    end

    prefix = toric3d_encoded_prefix(code, 1)
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    prefix_logical = first(toric3d_encoded_logical_x(prefix.code))
    frame.x[prefix.qubits[prefix_logical]] .= true
    checkpoint = TopoNoise._toric3d_apply_x_prefix_feedback!(
        frame, prefix, decoders.prefixes[1])
    @test !any(checkpoint.pre_syndrome)
    @test !any(checkpoint.recovery)
    @test !any(checkpoint.post_syndrome)
    propagate_toric3d_segment!(frame, encoder, clean,
        encoder.checkpoints[1] + 1, encoder.checkpoints[2])
    logical_outcome = toric3d_x_threshold_outcome(decoders.prefixes[end], frame.x)
    @test (logical_outcome.logical_1, logical_outcome.logical_2,
           logical_outcome.failure) == (true, false, true)
end

@testset "X strength edge cases and validation" begin
    one = Toric3DEncodedSlice(1, 2, 2)
    one_encoder = toric3d_encoded_encoder(one)
    one_decoders = build_toric3d_x_feedback_decoders(one; p=0.02)
    one_faults = sample_toric3d_faults(MersenneTwister(31), one_encoder;
        p_x=0.2, p_z=0.0)
    @test replay_toric3d_x_strength(one_encoder, one_faults, one_decoders;
        protocol=:full_prefix, record_history=true) ==
          replay_toric3d_x_strength(one_encoder, one_faults, one_decoders;
              protocol=:none, record_history=true)

    code = Toric3DEncodedSlice(2, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    zero_prior = build_toric3d_x_feedback_decoders(code; p=0.0)
    @test all(decoder.prior == 1e-12 for decoder in zero_prior.prefixes)
    prefix_decoder = zero_prior.prefixes[1]
    data_error = falses(size(prefix_decoder.matrix, 2))
    data_error[1] = true
    syndrome = toric3d_strength_x_syndrome(prefix_decoder.matrix, data_error)
    @test any(syndrome)
    @test toric3d_strength_x_syndrome(prefix_decoder.matrix,
        decode_toric3d_x(prefix_decoder, syndrome)) == syndrome

    clean = sample_toric3d_faults(MersenneTwister(37), encoder)
    @test_throws ArgumentError replay_toric3d_x_strength(
        encoder, clean, zero_prior; protocol=:bad)
    @test_throws ArgumentError replay_toric3d_x_strength(
        encoder, clean, zero_prior; protocol=:none, routing=:bad)
    wrong_encoder = toric3d_encoded_encoder(Toric3DEncodedSlice(1, 2, 2))
    wrong_faults = sample_toric3d_faults(MersenneTwister(41), wrong_encoder)
    @test_throws DimensionMismatch replay_toric3d_x_strength(
        wrong_encoder, wrong_faults, zero_prior; protocol=:none)
    wrong_decoders = build_toric3d_x_feedback_decoders(
        Toric3DEncodedSlice(1, 2, 2); p=0.0)
    @test_throws DimensionMismatch replay_toric3d_x_strength(
        encoder, clean, wrong_decoders; protocol=:full_prefix)
    z_fault = with_toric3d_fault(encoder, clean, 1,
        only(encoder.operations[1].qubits), :Z)
    @test_throws ArgumentError replay_toric3d_x_strength(
        encoder, z_fault, zero_prior; protocol=:full_prefix)
    @test_throws ArgumentError build_toric3d_x_feedback_decoders(code; p=-0.1)
end
