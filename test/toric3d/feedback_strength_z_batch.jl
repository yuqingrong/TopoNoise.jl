using Test
using Random
using TopoNoise

function toric3d_z_batch_scalar_record(batch, shot)
    transfer = Toric3DFaultRecord(batch.steps, batch.Ly, batch.Lz,
        falses(size(batch.transfer, 1), size(batch.transfer, 2)),
        BitMatrix(batch.transfer[:, :, shot]))
    return Toric3DDiagonalZFaultRecord(batch.seed[:, shot], transfer)
end

@testset "Batched Z sampling exactly matches scalar draw order" begin
    code = Toric3DEncodedSlice(2, 2, 3)
    encoder = toric3d_encoded_encoder(code)
    for shots in (1, 7, 16), p in (0.0, 0.27, 1.0)
        batch_rng = Xoshiro(101 + shots)
        scalar_rng = Xoshiro(101 + shots)
        batch = sample_toric3d_diagonal_z_fault_batch(batch_rng, encoder; p, shots)
        @test batch isa Toric3DDiagonalZFaultBatch
        @test (batch.steps, batch.Ly, batch.Lz) == (2, 2, 3)
        @test size(batch.seed) == (18, shots)
        @test size(batch.transfer) == (2, length(encoder.operations), shots)
        for shot in 1:shots
            scalar = sample_toric3d_diagonal_z_faults(scalar_rng, encoder; p)
            @test batch.seed[:, shot] == scalar.seed
            @test batch.transfer[:, :, shot] == scalar.transfer.z
            @test toric3d_z_batch_scalar_record(batch, shot) == scalar
        end
    end

    @test_throws ArgumentError sample_toric3d_diagonal_z_fault_batch(
        Xoshiro(1), encoder; p=-0.1, shots=2)
    @test_throws ArgumentError sample_toric3d_diagonal_z_fault_batch(
        Xoshiro(1), encoder; p=1.1, shots=2)
    @test_throws ArgumentError sample_toric3d_diagonal_z_fault_batch(
        Xoshiro(1), encoder; p=0.1, shots=0)
end

@testset "Batched Z replay matches scalar protocols column by column" begin
    code = Toric3DEncodedSlice(2, 2, 3)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    faults = sample_toric3d_diagonal_z_fault_batch(Xoshiro(211), encoder;
        p=0.19, shots=17)
    saved_seed = copy(faults.seed)
    saved_transfer = copy(faults.transfer)
    batch = toric3d_z_strength_batch(encoder, faults, decoders)
    @test keys(batch) == (:none, :slice_local, :full_prefix)

    for shot in 1:17
        scalar_faults = toric3d_z_batch_scalar_record(faults, shot)
        scalar = toric3d_z_strength_triplet(encoder, scalar_faults, decoders)
        for protocol in keys(batch)
            batch_result = getproperty(batch, protocol)
            scalar_result = getproperty(scalar, protocol)
            @test batch_result.logical_1[shot] == scalar_result.logical_1
            @test batch_result.logical_2[shot] == scalar_result.logical_2
            @test batch_result.failure[shot] == scalar_result.failure
        end
    end
    @test faults.seed == saved_seed
    @test faults.transfer == saved_transfer

    repeated = sample_toric3d_diagonal_z_fault_batch(Xoshiro(211), encoder;
        p=0.19, shots=17)
    @test faults == repeated
    @test toric3d_z_strength_batch(encoder, repeated, decoders) == batch

    zero = sample_toric3d_diagonal_z_fault_batch(Xoshiro(223), encoder;
        p=0.0, shots=7)
    zero_result = toric3d_z_strength_batch(encoder, zero, decoders)
    @test all(!any(result.failure) for result in values(zero_result))

    matrix = BitMatrix(Matrix(decoders.transfer.matrix) .!= 0)
    recovery = decode_toric3d_z_feedback_batch(decoders.transfer, matrix)
    @test TopoNoise._toric3d_validate_z_feedback_batch(
        decoders.transfer, matrix, recovery) == recovery
    corrupted = copy(recovery)
    corrupted[1, 1] ⊻= true
    @test_throws ErrorException TopoNoise._toric3d_validate_z_feedback_batch(
        decoders.transfer, matrix, corrupted)

    wrong = Toric3DEncodedSlice(1, 2, 3)
    wrong_encoder = toric3d_encoded_encoder(wrong)
    wrong_faults = sample_toric3d_diagonal_z_fault_batch(Xoshiro(227),
        wrong_encoder; p=0.1, shots=3)
    @test_throws DimensionMismatch toric3d_z_strength_batch(
        wrong_encoder, wrong_faults, decoders)
end

@testset "One-transfer batched local and prefix feedback coincide" begin
    code = Toric3DEncodedSlice(1, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    faults = sample_toric3d_diagonal_z_fault_batch(Xoshiro(229), encoder;
        p=0.23, shots=11)
    result = toric3d_z_strength_batch(encoder, faults, decoders)
    @test result.slice_local == result.full_prefix
end

@testset "Batched Z feedback measures the complete three-sheet window" begin
    code = Toric3DEncodedSlice(3, 2, 2)
    encoder = toric3d_encoded_encoder(code)
    decoders = build_toric3d_z_feedback_decoders(code)
    seed = falses(3code.lattice.Ly * code.lattice.Lz, 1)
    transfer = falses(2, length(encoder.operations), 1)
    transfer[2, 2, 1] = true
    faults = Toric3DDiagonalZFaultBatch(
        code.lattice.steps, code.lattice.Ly, code.lattice.Lz, seed, transfer)

    scalar = toric3d_z_strength_triplet(
        encoder, toric3d_z_batch_scalar_record(faults, 1), decoders)
    batch = toric3d_z_strength_batch(encoder, faults, decoders)
    for protocol in keys(batch)
        batch_result = getproperty(batch, protocol)
        scalar_result = getproperty(scalar, protocol)
        @test only(batch_result.logical_1) == scalar_result.logical_1
        @test only(batch_result.logical_2) == scalar_result.logical_2
        @test only(batch_result.failure) == scalar_result.failure
    end
end
