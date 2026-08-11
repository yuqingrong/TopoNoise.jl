using TopoNoise
using Random
using Test

@testset "open code-capacity data edges" begin
    model = OpenCodeCapacityModel(4, 5)
    zero = sample_data_edge_errors(MersenneTwister(1), model; error_rate=0.0)
    one = sample_data_edge_errors(MersenneTwister(1), model; error_rate=1.0)
    @test size(zero.horizontal) == (4, 4)
    @test size(zero.vertical) == (3, 5)
    @test !any(zero.horizontal) && !any(zero.vertical)
    @test all(one.horizontal) && all(one.vertical)
    @test !any(code_capacity_syndrome(zero))
end

@testset "one edge has the plaquette-incidence syndrome" begin
    errors = DataEdgeErrors(falses(4, 3), falses(3, 4))
    errors.horizontal[2, 2] = true
    syndrome = code_capacity_syndrome(errors)
    @test syndrome[1, 2] && syndrome[2, 2]
    @test count(syndrome) == 2
end

@testset "data-edge inputs and logical cuts are validated" begin
    @test_throws ArgumentError OpenCodeCapacityModel(1, 5)
    @test_throws ArgumentError OpenCodeCapacityModel(4, 1)
    @test_throws ArgumentError DataEdgeErrors(falses(4, 4), falses(2, 5))

    model = OpenCodeCapacityModel(4, 5)
    north_south = logical_cut(model; sector=:north_south)
    east_west = logical_cut(model; sector=:east_west)
    @test north_south.horizontal == BitMatrix([false false false false;
                                               true true true true;
                                               false false false false;
                                               false false false false])
    @test !any(north_south.vertical)
    @test east_west.vertical == BitMatrix([false false true false false;
                                           false false true false false;
                                           false false true false false])
    @test !any(east_west.horizontal)
    @test_throws ArgumentError logical_cut(model; sector=:diagonal)
end

function single_edge_error_configurations(model)
    errors = DataEdgeErrors[]
    for r in 1:model.rows, c in 1:(model.cols - 1)
        horizontal = falses(model.rows, model.cols - 1)
        horizontal[r, c] = true
        push!(errors, DataEdgeErrors(horizontal, falses(model.rows - 1, model.cols)))
    end
    for r in 1:(model.rows - 1), c in 1:model.cols
        vertical = falses(model.rows - 1, model.cols)
        vertical[r, c] = true
        push!(errors, DataEdgeErrors(falses(model.rows, model.cols - 1), vertical))
    end
    errors
end

@testset "syndrome-only decode" begin
    model = OpenCodeCapacityModel(4, 4)
    for sector in (:north_south, :east_west)
        for errors in single_edge_error_configurations(model)
            syndrome = code_capacity_syndrome(errors)
            correction = decode_syndrome(model, syndrome; sector=sector)
            residual = residual_errors(errors, correction)
            @test !any(code_capacity_syndrome(residual))
            @test !logical_failure(errors, correction; sector=sector)
        end
    end
end

@testset "logical-cut geometry" begin
    model = OpenCodeCapacityModel(4, 4)
    none = Correction(4, 4)
    north_south = DataEdgeErrors(falses(4, 3), falses(3, 4))
    north_south.horizontal[:, 2] .= true
    @test !any(code_capacity_syndrome(north_south))
    @test logical_failure(north_south, none; sector=:north_south)
    @test !logical_failure(north_south, none; sector=:east_west)

    east_west = DataEdgeErrors(falses(4, 3), falses(3, 4))
    east_west.vertical[2, :] .= true
    @test !any(code_capacity_syndrome(east_west))
    @test logical_failure(east_west, none; sector=:east_west)
    @test !logical_failure(east_west, none; sector=:north_south)
end
