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
    @test north_south.vertical == BitMatrix([false false false false false;
                                             true true true true true;
                                             false false false false false])
    @test !any(north_south.horizontal)
    @test east_west.horizontal == BitMatrix([false true false false;
                                             false true false false;
                                             false true false false;
                                             false true false false])
    @test !any(east_west.vertical)
    @test_throws ArgumentError logical_cut(model; sector=:diagonal)
end
