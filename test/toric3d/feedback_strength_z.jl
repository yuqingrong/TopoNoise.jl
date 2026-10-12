using Test
using Random
using TopoNoise

toric3d_strength_z_syndrome(matrix, bits) =
    BitVector(isodd.(matrix * Vector{UInt8}(bits)))

function toric3d_strength_z_validate_final(result, decoder)
    @test toric3d_strength_z_syndrome(decoder.matrix, result.recovery) == result.syndrome
    @test !any(toric3d_strength_z_syndrome(decoder.matrix, result.residual))
    expected = toric3d_z_threshold_outcome(decoder, result.frame.z)
    @test result.syndrome == expected.syndrome
    @test result.recovery == expected.recovery
    @test result.residual == expected.residual
    @test (result.logical_1, result.logical_2, result.failure) ==
          (expected.logical_1, expected.logical_2, expected.failure)
end

@testset "Diagonal Z fault records contain every eligible location in order" begin
    code = Toric3DEncodedSlice(2, 2, 3)
    encoder = toric3d_encoded_encoder(code)
    cells = code.lattice.Ly * code.lattice.Lz
    locations = sum(length(operation.qubits) for operation in encoder.operations)

    first_record = sample_toric3d_diagonal_z_faults(MersenneTwister(17), encoder; p=0.23)
    second_record = sample_toric3d_diagonal_z_faults(MersenneTwister(17), encoder; p=0.23)
    @test first_record isa Toric3DDiagonalZFaultRecord
    @test first_record == second_record
    @test length(first_record.seed) == 3cells
    @test size(first_record.transfer.z) == (2, length(encoder.operations))
    @test !any(first_record.transfer.x)
    @test count(first_record.seed) + count(first_record.transfer.z) <= 3cells + locations

    zero = sample_toric3d_diagonal_z_faults(MersenneTwister(19), encoder; p=0.0)
    one = sample_toric3d_diagonal_z_faults(MersenneTwister(19), encoder; p=1.0)
    @test !any(zero.seed)
    @test !any(zero.transfer.x)
    @test !any(zero.transfer.z)
    @test all(one.seed)
    @test !any(one.transfer.x)
    @test count(one.transfer.z) == locations
    @test all(one.transfer.z[slot, step]
              for (step, operation) in enumerate(encoder.operations)
              for slot in eachindex(operation.qubits))

    @test_throws ArgumentError sample_toric3d_diagonal_z_faults(
        MersenneTwister(1), encoder; p=-0.1)
    @test_throws ArgumentError sample_toric3d_diagonal_z_faults(
        MersenneTwister(1), encoder; p=1.1)
    @test_throws DimensionMismatch Toric3DDiagonalZFaultRecord(
        falses(3cells - 1), zero.transfer)
    x_transfer = with_toric3d_fault(encoder, zero.transfer, 1,
        only(encoder.operations[1].qubits), :X)
    @test_throws ArgumentError Toric3DDiagonalZFaultRecord(falses(3cells), x_transfer)
end

@testset "Local and full-prefix Z feedback clear every completed prefix" begin
    code = Toric3DEncodedSlice(3, 2, 3)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    record = sample_toric3d_diagonal_z_faults(MersenneTwister(23), encoder; p=0.19)
    saved_seed = copy(record.seed)
    saved_transfer_x = copy(record.transfer.x)
    saved_transfer_z = copy(record.transfer.z)
    triplet = toric3d_z_strength_triplet(encoder, record, decoders;
        record_history=true)

    @test keys(triplet) == (:none, :slice_local, :full_prefix)
    @test isempty(triplet.none.history)
    @test length(triplet.slice_local.history) == 3
    @test length(triplet.full_prefix.history) == 3
    @test triplet.slice_local.history[1].kind === :seed
    @test triplet.full_prefix.history[1].kind === :seed
    @test triplet.slice_local.history[1].syndrome ==
          triplet.full_prefix.history[1].syndrome
    @test triplet.slice_local.history[1].recovery ==
          triplet.full_prefix.history[1].recovery
    @test triplet.slice_local.history[1].qubits ==
          triplet.full_prefix.history[1].qubits

    for checkpoint in triplet.slice_local.history
        @test !any(checkpoint.post_syndrome)
        @test checkpoint.support == checkpoint.qubits[findall(checkpoint.recovery)]
        if checkpoint.checkpoint > 0
            expected = toric3d_z_feedback_window(code, checkpoint.checkpoint).qubits
            @test checkpoint.qubits == expected
            @test all(data_qubit_coordinate(code, q)[2] in
                      (checkpoint.checkpoint - 1:checkpoint.checkpoint)
                      for q in checkpoint.support)
        end
    end
    for checkpoint in triplet.full_prefix.history
        @test !any(checkpoint.post_syndrome)
        @test checkpoint.support == checkpoint.qubits[findall(checkpoint.recovery)]
        if checkpoint.checkpoint > 0
            expected = toric3d_encoded_prefix(code, checkpoint.checkpoint).qubits
            @test checkpoint.qubits == expected
            @test all(data_qubit_coordinate(code, q)[2] <= checkpoint.checkpoint
                      for q in checkpoint.support)
        end
    end

    for result in values(triplet)
        @test !any(result.frame.x)
        toric3d_strength_z_validate_final(result, decoders.prefixes[end])
    end
    @test record.seed == saved_seed
    @test record.transfer.x == saved_transfer_x
    @test record.transfer.z == saved_transfer_z
end

@testset "Every small Z fault has a valid three-protocol outcome" begin
    code = Toric3DEncodedSlice(2, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    clean = sample_toric3d_diagonal_z_faults(MersenneTwister(29), encoder; p=0.0)
    records = Toric3DDiagonalZFaultRecord[]
    for seed_qubit in eachindex(clean.seed)
        seed = copy(clean.seed)
        seed[seed_qubit] = true
        push!(records, Toric3DDiagonalZFaultRecord(seed, clean.transfer))
    end
    for (gate, operation) in enumerate(encoder.operations), qubit in operation.qubits
        transfer = with_toric3d_fault(encoder, clean.transfer, gate, qubit, :Z)
        push!(records, Toric3DDiagonalZFaultRecord(clean.seed, transfer))
    end
    @test length(records) == 3 * 2 * 2 + 13 * 2 * 2 * 2

    for record in records
        triplet = toric3d_z_strength_triplet(encoder, record, decoders;
            record_history=true)
        for result in values(triplet)
            toric3d_strength_z_validate_final(result, decoders.prefixes[end])
        end
        @test all(!any(checkpoint.post_syndrome)
                  for checkpoint in triplet.slice_local.history)
        @test all(!any(checkpoint.post_syndrome)
                  for checkpoint in triplet.full_prefix.history)
    end
end

@testset "Three-sheet Z window contains all transfer-fault endpoints" begin
    code = Toric3DEncodedSlice(1, 3, 3)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    clean = sample_toric3d_diagonal_z_faults(MersenneTwister(31), encoder; p=0.0)
    window = toric3d_z_feedback_window(code, 1)
    cells = 9
    odd_two_sheet = 0
    checked = 0
    for (gate, operation) in enumerate(encoder.operations), qubit in operation.qubits
        transfer = with_toric3d_fault(encoder, clean.transfer, gate, qubit, :Z)
        frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
        propagate_toric3d_segment!(frame, encoder, transfer, 1, length(encoder.operations))
        syndrome = toric3d_z_window_syndrome(window, frame.z)
        odd_two_sheet += isodd(count(@view syndrome[cells + 1:end]))
        @test iseven(count(syndrome))
        recovery = decode_toric3d_z_feedback(decoders.transfer, syndrome)
        @test toric3d_strength_z_syndrome(decoders.transfer.matrix, recovery) == syndrome
        checked += 1
    end
    @test checked == 117
    @test odd_two_sheet == 54
end

@testset "Random Z histories and edge cases preserve logical information" begin
    code = Toric3DEncodedSlice(3, 2, 3)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    for seed in 1:8
        record = sample_toric3d_diagonal_z_faults(MersenneTwister(seed), encoder; p=0.22)
        triplet = toric3d_z_strength_triplet(encoder, record, decoders;
            record_history=true)
        @test all(!any(checkpoint.post_syndrome)
                  for checkpoint in triplet.slice_local.history)
        @test all(!any(checkpoint.post_syndrome)
                  for checkpoint in triplet.full_prefix.history)
    end

    one = Toric3DEncodedSlice(1, 2, 2)
    one_encoder = toric3d_encoded_encoder(one)
    one_decoders = build_toric3d_z_feedback_decoders(one)
    one_faults = sample_toric3d_diagonal_z_faults(MersenneTwister(37), one_encoder; p=0.2)
    one_triplet = toric3d_z_strength_triplet(one_encoder, one_faults, one_decoders;
        record_history=true)
    @test one_triplet.slice_local == one_triplet.full_prefix

    clean = sample_toric3d_diagonal_z_faults(MersenneTwister(41), encoder; p=0.0)
    logical_seed = copy(clean.seed)
    logical_seed[first(toric3d_encoded_logical_z(code))] .= true
    logical_record = Toric3DDiagonalZFaultRecord(logical_seed, clean.transfer)
    logical_result = replay_toric3d_z_strength(encoder, logical_record, decoders;
        protocol=:slice_local, record_history=true)
    @test !any(logical_result.history[1].syndrome)
    @test !any(logical_result.history[1].recovery)
    @test (logical_result.logical_1, logical_result.logical_2,
           logical_result.failure) == (true, false, true)

    final_gate = encoder.checkpoints[end - 1] + 1
    final_transfer = with_toric3d_fault(encoder, clean.transfer, final_gate,
        only(encoder.operations[final_gate].qubits), :Z)
    final_record = Toric3DDiagonalZFaultRecord(clean.seed, final_transfer)
    final_triplet = toric3d_z_strength_triplet(encoder, final_record, decoders)
    @test final_triplet.none == final_triplet.slice_local == final_triplet.full_prefix
    @test any(final_triplet.none.syndrome)

    bad_syndrome = falses(size(decoders.transfer.matrix, 1))
    bad_syndrome[1:2] .= true
    @test_throws ErrorException TopoNoise._toric3d_validate_z_feedback_recovery(
        decoders.transfer, bad_syndrome, falses(size(decoders.transfer.matrix, 2)))
    @test_throws ArgumentError replay_toric3d_z_strength(
        encoder, clean, decoders; protocol=:bad)
    wrong_code = Toric3DEncodedSlice(2, 2, 3)
    wrong_encoder = toric3d_encoded_encoder(wrong_code)
    wrong_record = sample_toric3d_diagonal_z_faults(MersenneTwister(43), wrong_encoder; p=0.0)
    @test_throws DimensionMismatch replay_toric3d_z_strength(
        wrong_encoder, wrong_record, decoders; protocol=:none)
    wrong_decoders = build_toric3d_z_feedback_decoders(wrong_code)
    @test_throws DimensionMismatch replay_toric3d_z_strength(
        encoder, clean, wrong_decoders; protocol=:none)
end
