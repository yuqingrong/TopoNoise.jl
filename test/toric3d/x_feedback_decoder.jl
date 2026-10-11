using Test
using TopoNoise

function toric3d_x_feedback_test_syndrome(matrix, error)
    return BitVector(isodd(count((@view matrix[row, :]) .& error))
        for row in axes(matrix, 1))
end

@testset "Encoded diagonal X recovery is syndrome-valid and scores logicals" begin
    code = Toric3DEncodedSlice(2, 2, 3)
    n = data_qubit_count(code)
    matrix = BitMatrix(plaquette_check_matrix(code))
    settings = Toric3DDecoderSettings()
    decoder = build_toric3d_x_threshold_decoder(code; p=0.0, settings)

    @test decoder.code === code
    @test decoder.p == 0.0
    @test decoder.prior == 1e-12
    @test decoder.settings === settings
    @test decoder.matrix == matrix
    @test size(decoder.logicals) == (2, n)
    @test !isempty(decoder.version)

    for q in 1:n
        error = falses(n)
        error[q] = true
        syndrome = toric3d_x_feedback_test_syndrome(matrix, error)
        recovery = decode_toric3d_x(decoder, syndrome)
        @test length(recovery) == n
        @test toric3d_x_feedback_test_syndrome(matrix, recovery) == syndrome
    end

    logicals = toric3d_encoded_logical_x(code)
    for logical_index in 1:2
        frame = falses(n)
        frame[logicals[logical_index]] .= true
        outcome = toric3d_x_threshold_outcome(decoder, frame)
        expected = logical_index == 1 ? (true, false, true) : (false, true, true)
        @test (outcome.logical_1, outcome.logical_2, outcome.failure) == expected
        @test !any(outcome.syndrome)
        @test !any(outcome.recovery)
        @test outcome.residual == frame
    end

    stabilizer = falses(n)
    stabilizer[first(star_checks(code))] .= true
    stable = toric3d_x_threshold_outcome(decoder, stabilizer)
    @test !stable.logical_1
    @test !stable.logical_2
    @test !stable.failure
    @test !any(stable.syndrome)

    nonzero = falses(n)
    nonzero[1] = true
    nonzero_syndrome = toric3d_x_feedback_test_syndrome(matrix, nonzero)
    @test any(nonzero_syndrome)
    @test toric3d_x_feedback_test_syndrome(matrix,
        decode_toric3d_x(decoder, nonzero_syndrome)) == nonzero_syndrome

    @test_throws DimensionMismatch decode_toric3d_x(decoder,
        falses(size(matrix, 1) - 1))
    @test_throws DimensionMismatch TopoNoise._toric3d_validate_x_recovery(
        decoder, nonzero_syndrome, falses(n - 1))
    @test_throws ErrorException TopoNoise._toric3d_validate_x_recovery(
        decoder, nonzero_syndrome, falses(n))
    @test_throws DimensionMismatch toric3d_x_threshold_outcome(decoder, falses(n - 1))
    @test_throws ArgumentError build_toric3d_x_threshold_decoder(code; p=-0.1)
    @test_throws ArgumentError build_toric3d_x_threshold_decoder(code; p=1.1)

    numpy = TopoNoise.PythonCall.pyimport("numpy")
    @test_throws ErrorException TopoNoise._toric3d_py_bits(
        numpy, numpy.asarray(UInt8[0, 1]), 3)
    @test_throws ErrorException TopoNoise._toric3d_py_bits(
        numpy, numpy.asarray(UInt8[0, 2, 1]), 3)
end
