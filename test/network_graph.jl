using ITensorNetworks: ITensorNetworks, ITensorNetwork

function canonical_edge(edge)
    first_vertex, second_vertex = Tuple(edge)
    return isless(first_vertex, second_vertex) ?
           (first_vertex, second_vertex) : (second_vertex, first_vertex)
end

@testset "Toric-code PEPS graph" begin
    @test isdefined(TopoNoise, :peps_graph)

    if isdefined(TopoNoise, :peps_graph)
        for (rows, cols) in ((1, 1), (1, 2), (2, 1), (2, 2), (3, 3))
            sized_peps = toric_code_peps(rows, cols)
            sized_network = peps_graph(sized_peps)
            expected_vertices =
                Set((row, col) for row in 1:rows for col in 1:cols)
            expected_edges = rows * (cols - 1) + (rows - 1) * cols
            @test Set(keys(sized_network)) == expected_vertices
            @test length(collect(ITensorNetworks.edges(sized_network))) == expected_edges
            @test all(hasinds(sized_network[(row, col)],
                              physicalinds(sized_peps, row, col)...)
                      for row in 1:rows for col in 1:cols)
        end

        peps = toric_code_peps(3, 3)
        network = TopoNoise.peps_graph(peps)

        @test network isa ITensorNetwork{Tuple{Int,Int}}
        @test Set(keys(network)) == Set((row, col) for row in 1:3 for col in 1:3)

        expected_edges = Set{Tuple{Tuple{Int,Int},Tuple{Int,Int}}}()
        for row in 1:3, col in 1:2
            push!(expected_edges, ((row, col), (row, col + 1)))
        end
        for row in 1:2, col in 1:3
            push!(expected_edges, ((row, col), (row + 1, col)))
        end
        actual_edges = Set(canonical_edge(edge) for edge in ITensorNetworks.edges(network))
        @test actual_edges == expected_edges

        for row in 1:3, col in 1:3
            @test network[(row, col)] == peps[row, col]
            @test hasinds(network[(row, col)], physicalinds(peps, row, col)...)
        end

        malformed_shape = ToricCodePEPS(
            copy(peps.tensors), copy(peps.physical_indices[1:2, :]))
        @test_throws DimensionMismatch TopoNoise.peps_graph(malformed_shape)

        duplicated_physical = copy(peps.physical_indices)
        east, north, west, south = duplicated_physical[1, 1]
        duplicated_physical[1, 1] = (east, east, west, south)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(copy(peps.tensors), duplicated_physical))

        misordered_physical = copy(peps.physical_indices)
        misordered_physical[1, 1] = (north, east, west, south)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(copy(peps.tensors), misordered_physical))

        wrong_site_index = Index(2, "Site,east,r=99,c=99")
        wrong_site_tensors = copy(peps.tensors)
        wrong_site_tensors[1, 1] = replaceind(
            wrong_site_tensors[1, 1], east, wrong_site_index)
        wrong_site_physical = copy(peps.physical_indices)
        wrong_site_physical[1, 1] =
            (wrong_site_index, north, west, south)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(wrong_site_tensors, wrong_site_physical))

        ambiguous_direction = Index(2, "Site,east,north,r=1,c=1")
        ambiguous_tensors = copy(peps.tensors)
        ambiguous_tensors[1, 1] = replaceind(
            ambiguous_tensors[1, 1], east, ambiguous_direction)
        ambiguous_physical = copy(peps.physical_indices)
        ambiguous_physical[1, 1] =
            (ambiguous_direction, north, west, south)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(ambiguous_tensors, ambiguous_physical))

        bad_physical_dimension = Index(3, "Site,east,r=1,c=1")
        one_site = toric_code_peps(1, 1)
        _, one_north, one_west, one_south = physicalinds(one_site, 1, 1)
        bad_physical_tensor = ITensor(
            Float64, bad_physical_dimension, one_north, one_west, one_south)
        bad_physical_peps = ToricCodePEPS(
            reshape([bad_physical_tensor], 1, 1),
            reshape([(bad_physical_dimension, one_north, one_west, one_south)], 1, 1))
        @test_throws ArgumentError peps_graph(bad_physical_peps)
        @test_throws ArgumentError sequential_circuit_graph(bad_physical_peps)

        two_sites = toric_code_peps(1, 2)
        bad_virtual_dimension = Index(3, "Link,h,r=1,c=1")
        bad_virtual_tensors = Matrix{ITensor}(undef, 1, 2)
        for col in 1:2
            bad_virtual_tensors[1, col] = ITensor(
                Float64, physicalinds(two_sites, 1, col)...,
                bad_virtual_dimension)
        end
        bad_virtual_peps = ToricCodePEPS(
            bad_virtual_tensors, copy(two_sites.physical_indices))
        @test_throws ArgumentError peps_graph(bad_virtual_peps)
        @test_throws ArgumentError sequential_circuit_graph(bad_virtual_peps)

        extra_bond = Index(2, "Link,unexpected")
        diagonal_tensors = copy(peps.tensors)
        diagonal_tensors[1, 1] *= onehot(extra_bond => 1)
        diagonal_tensors[2, 2] *= onehot(extra_bond => 1)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(diagonal_tensors, copy(peps.physical_indices)))

        hyperedge = Index(2, "Link,hyperedge")
        hyperedge_tensors = copy(peps.tensors)
        for site in ((1, 1), (1, 2), (2, 1))
            hyperedge_tensors[site...] *= onehot(hyperedge => 1)
        end
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(hyperedge_tensors, copy(peps.physical_indices)))

        empty_tensors = Matrix{ITensor}(undef, 0, 0)
        empty_physical = Matrix{NTuple{4,Index}}(undef, 0, 0)
        @test_throws ArgumentError TopoNoise.peps_graph(
            ToricCodePEPS(empty_tensors, empty_physical))
    end
end
