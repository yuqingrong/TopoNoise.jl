using Random
using Test

function _planar_zero_errors(model)
    IsometricPlanarVirtualErrors(
        falses(model.rows, model.cols - 1), falses(model.rows - 1, model.cols))
end

function _single_planar_error(model, direction::Symbol, row::Int, col::Int)
    errors = _planar_zero_errors(model)
    getproperty(errors, direction)[row, col] = true
    return errors
end

@testset "Pseudo-toric API is retired" begin
    @test !isdefined(TopoNoise, :IsometricToricCodeModel)
    @test !isdefined(TopoNoise, :scan_isometric_toric_capacity)
end

@testset "Isometric planar-code geometry and detectors" begin
    model = IsometricPlanarCodeModel(3)
    @test model.distance == 3
    @test (model.rows, model.cols) == (3, 3)

    errors = IsometricPlanarVirtualErrors(falses(3, 2), falses(2, 3))
    syndrome = isometric_planar_syndrome(errors)
    @test size(syndrome.plaquettes) == (2, 2)
    @test length(syndrome.north) == 2
    @test length(syndrome.south) == 2
    @test !any(syndrome.plaquettes)
    @test !any(syndrome.north)
    @test !any(syndrome.south)
end

@testset "Isometric planar-code Yao detector reference" begin
    for distance in (2, 3)
        model = IsometricPlanarCodeModel(distance)
        error_records = [_planar_zero_errors(model)]
        for row in 1:model.rows, col in 1:(model.cols - 1)
            push!(error_records, _single_planar_error(model, :horizontal, row, col))
        end
        for row in 1:(model.rows - 1), col in 1:model.cols
            push!(error_records, _single_planar_error(model, :vertical, row, col))
        end
        append!(error_records, [
            sample_isometric_planar_virtual_errors(
                MersenneTwister(1000 + distance + sample), model; error_rate=0.2)
            for sample in 1:3])

        for (index, errors) in enumerate(error_records)
            trajectory = sample_isometric_planar_yao_trajectory(
                MersenneTwister(2000 + 100 * distance + index), model, errors)
            yao_syndrome = isometric_planar_syndrome(
                isometric_planar_bond_mismatches(trajectory))
            direct_syndrome = isometric_planar_syndrome(errors)
            @test yao_syndrome.plaquettes == direct_syndrome.plaquettes
            @test yao_syndrome.north == direct_syndrome.north
            @test yao_syndrome.south == direct_syndrome.south
            @test TopoNoise._isometric_planar_yao_logical_frame(trajectory) ==
                  TopoNoise._isometric_planar_logical_frame(errors)
        end
    end
end

