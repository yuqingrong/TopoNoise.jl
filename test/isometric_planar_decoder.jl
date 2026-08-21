using Random
using Test

@testset "Isometric planar-code encoder and parity-check layout" begin
    model = IsometricPlanarCodeModel(3)
    encoder = isometric_planar_encoder(model)
    @test encoder.wire_count == 36
    @test size(encoder.output_wires) == (3, 3, 4)
    @test sort(vec(encoder.output_wires)) == collect(1:36)

    checks = isometric_planar_check_layer(encoder)
    @test length(checks.detector_supports) == 8
    @test all(length(support) == 8 for support in checks.detector_supports[1:4])
    @test all(length(support) == 2 for support in checks.detector_supports[5:8])
    @test length(checks.logical_z_support) == 4
end

@testset "Isometric planar-code topology and decoder" begin
    for distance in 2:8
        model = IsometricPlanarCodeModel(distance)
        @test TopoNoise._isometric_planar_minimum_odd_logical_weight(model) == distance
    end

    model = IsometricPlanarCodeModel(4)
    logical = _planar_zero_errors(model)
    logical.vertical[1, :] .= true
    logical_syndrome = isometric_planar_syndrome(logical)
    @test !any(logical_syndrome.plaquettes)
    @test !any(logical_syndrome.north)
    @test !any(logical_syndrome.south)
    @test isometric_planar_logical_failure(logical, Correction(4, 4))

    corner = _planar_zero_errors(model)
    corner.horizontal[1, 1] = true
    corner.vertical[1, 1] = true
    corner_syndrome = isometric_planar_syndrome(corner)
    @test corner_syndrome.north == BitVector([true, false, false])
    @test !any(corner_syndrome.plaquettes)

    for distance in (3, 4)
        model = IsometricPlanarCodeModel(distance)
        for row in 1:model.rows, col in 1:(model.cols - 1)
            errors = _single_planar_error(model, :horizontal, row, col)
            correction = decode_isometric_planar_syndrome(
                model, isometric_planar_syndrome(errors))
            @test !isometric_planar_logical_failure(errors, correction)
        end
        for row in 1:(model.rows - 1), col in 1:model.cols
            errors = _single_planar_error(model, :vertical, row, col)
            correction = decode_isometric_planar_syndrome(
                model, isometric_planar_syndrome(errors))
            @test !isometric_planar_logical_failure(errors, correction)
        end
    end

    diagnostic = IsometricPlanarCodeModel(2)
    west = _single_planar_error(diagnostic, :vertical, 1, 1)
    correction = decode_isometric_planar_syndrome(
        diagnostic, isometric_planar_syndrome(west))
    @test correction.vertical[1, 1]
end

@testset "Isometric planar-code distance-five low-weight decoding" begin
    model = IsometricPlanarCodeModel(5)
    supports = Tuple{Symbol,Int,Int}[]
    append!(supports, [(:horizontal, row, col)
                       for row in 1:model.rows, col in 1:(model.cols - 1)])
    append!(supports, [(:vertical, row, col)
                       for row in 1:(model.rows - 1), col in 1:model.cols])

    records = IsometricPlanarVirtualErrors[]
    for support in supports
        errors = _planar_zero_errors(model)
        getproperty(errors, support[1])[support[2], support[3]] = true
        push!(records, errors)
    end
    for first in 1:(length(supports) - 1), second in (first + 1):length(supports)
        errors = _planar_zero_errors(model)
        for support in (supports[first], supports[second])
            getproperty(errors, support[1])[support[2], support[3]] = true
        end
        push!(records, errors)
    end

    detector_batch = falses(length(records), 24)
    scalar_frames = falses(length(records))
    for (index, errors) in enumerate(records)
        syndrome = isometric_planar_syndrome(errors)
        detector_batch[index, :] .= TopoNoise._isometric_planar_syndrome_vector(
            model, syndrome) .!= 0
        correction = decode_isometric_planar_syndrome(model, syndrome)
        @test !isometric_planar_logical_failure(errors, correction)
        scalar_frames[index] = xor(correction.vertical[:, 1]...)
    end
    @test TopoNoise._isometric_planar_batch_predictions(model, detector_batch) == scalar_frames
end
