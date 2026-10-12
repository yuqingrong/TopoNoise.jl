using Test
using TopoNoise

function toric3d_feedback_test_syndrome(matrix, frame)
    return BitVector(isodd.(matrix * Vector{UInt8}(frame)))
end

@testset "Z feedback decoder family validates local recoveries" begin
    code = Toric3DEncodedSlice(3, 2, 3)
    family = build_toric3d_z_feedback_decoders(code)
    @test family isa Toric3DZFeedbackDecoderFamily
    @test family.code === code
    @test family.seed isa Toric3DZFeedbackDecoder
    @test family.transfer isa Toric3DZFeedbackDecoder
    @test family.seed.window.r == 0
    @test family.transfer.window.r == 1
    @test length(family.prefixes) == 3
    @test [decoder.code.lattice.steps for decoder in family.prefixes] == [1, 2, 3]
    @test !isempty(family.seed.version)
    @test family.seed.version == family.transfer.version

    first_window = toric3d_z_feedback_window(code, 1)
    second_window = toric3d_z_feedback_window(code, 2)
    @test first_window.matrix == second_window.matrix
    @test BitMatrix(Matrix(family.transfer.matrix) .!= 0) == first_window.matrix

    for decoder in (family.seed, family.transfer)
        matrix = BitMatrix(Matrix(decoder.matrix) .!= 0)
        for column in axes(matrix, 2)
            syndrome = matrix[:, column]
            recovery = decode_toric3d_z_feedback(decoder, syndrome)
            @test length(recovery) == size(matrix, 2)
            @test toric3d_feedback_test_syndrome(decoder.matrix, recovery) == syndrome
        end

        recoveries = decode_toric3d_z_feedback_batch(decoder, matrix)
        @test size(recoveries) == (size(matrix, 2), size(matrix, 2))
        @test BitMatrix(isodd.(decoder.matrix * Matrix{UInt8}(recoveries))) == matrix
        @test size(decode_toric3d_z_feedback_batch(
            decoder, falses(size(matrix, 1), 0))) == (size(matrix, 2), 0)

        @test_throws DimensionMismatch decode_toric3d_z_feedback(
            decoder, falses(size(matrix, 1) - 1))
        @test_throws DimensionMismatch decode_toric3d_z_feedback_batch(
            decoder, falses(size(matrix, 1) - 1, 1))
        odd = falses(size(matrix, 1))
        odd[1] = true
        @test_throws ArgumentError decode_toric3d_z_feedback(decoder, odd)
        @test_throws ArgumentError decode_toric3d_z_feedback_batch(
            decoder, reshape(odd, :, 1))
    end

    for dims in ((1, 2, 2), (2, 2, 3), (2, 3, 2))
        transverse = build_toric3d_z_feedback_decoders(Toric3DEncodedSlice(dims...))
        @test all(==(2), vec(sum(transverse.seed.matrix; dims=1)))
        @test all(==(2), vec(sum(transverse.transfer.matrix; dims=1)))
    end
end

@testset "Single-shot final Z outcome preserves logical class" begin
    code = Toric3DEncodedSlice(2, 2, 3)
    decoder = build_toric3d_z_threshold_decoder(code)
    n = data_qubit_count(code)

    clean = toric3d_z_threshold_outcome(decoder, falses(n))
    @test keys(clean) == (:syndrome, :recovery, :residual,
                          :logical_1, :logical_2, :failure)
    @test !any(clean.syndrome)
    @test !any(clean.recovery)
    @test !any(clean.residual)
    @test (clean.logical_1, clean.logical_2, clean.failure) == (false, false, false)

    stabilizer = falses(n)
    stabilizer[first(plaquette_checks(code))] .= true
    stable = toric3d_z_threshold_outcome(decoder, stabilizer)
    @test !any(stable.syndrome)
    @test stable.residual == stabilizer
    @test (stable.logical_1, stable.logical_2, stable.failure) == (false, false, false)

    for (logical_index, support) in enumerate(toric3d_encoded_logical_z(code))
        frame = falses(n)
        frame[support] .= true
        outcome = toric3d_z_threshold_outcome(decoder, frame)
        expected = logical_index == 1 ? (true, false, true) : (false, true, true)
        @test !any(outcome.syndrome)
        @test !any(outcome.recovery)
        @test outcome.residual == frame
        @test (outcome.logical_1, outcome.logical_2, outcome.failure) == expected
    end

    single = falses(n)
    single[1] = true
    outcome = toric3d_z_threshold_outcome(decoder, single)
    @test toric3d_feedback_test_syndrome(decoder.matrix, outcome.recovery) ==
          outcome.syndrome
    corrupted = copy(outcome.recovery)
    corrupted[1] ⊻= true
    @test_throws ErrorException TopoNoise._toric3d_validate_z_recovery(
        decoder, outcome.syndrome, corrupted)
    @test_throws DimensionMismatch toric3d_z_threshold_outcome(decoder, falses(n - 1))
end
