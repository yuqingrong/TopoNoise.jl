using Test
using TopoNoise

@testset "Encoded diagonal-slice public geometry" begin
    @test isdefined(TopoNoise, :Toric3DEncodedSlice)
    @test isdefined(TopoNoise, :toric3d_encoded_logical_x)
    @test isdefined(TopoNoise, :verify_toric3d_encoded_encoder)
end

function toric3d_encoded_test_rank(matrix)
    work = BitMatrix(matrix)
    pivot = 0
    for column in axes(work, 2)
        source = findfirst(row -> work[row, column], pivot+1:size(work, 1))
        isnothing(source) && continue
        source += pivot
        pivot += 1
        work[pivot, :], work[source, :] = copy(work[source, :]), copy(work[pivot, :])
        for row in pivot+1:size(work, 1)
            work[row, column] && (work[row, :] .⊻= work[pivot, :])
        end
    end
    return pivot
end

# Shortest closed syndrome-graph path with either logical winding odd. The
# four-sheet lift retains logical parity, including parallel length-two edges.
function toric3d_encoded_test_z_distance(code)
    hx = star_check_matrix(code)
    logicals = TopoNoise.toric3d_encoded_logical_x(code)
    graph = [Tuple{Int,Int}[] for _ in axes(hx, 1)]
    for q in axes(hx, 2)
        a, b = findall(hx[:, q])
        winding = Int(q in logicals[1]) + 2Int(q in logicals[2])
        push!(graph[a], (b, winding))
        push!(graph[b], (a, winding))
    end
    shortest = typemax(Int)
    for start in eachindex(graph)
        distances = fill(-1, length(graph), 4)
        distances[start, 1] = 0
        queue = [(start, 0)]
        cursor = 1
        while cursor <= length(queue)
            node, winding = queue[cursor]
            cursor += 1
            for (neighbor, edge_winding) in graph[node]
                next_winding = xor(winding, edge_winding)
                distances[neighbor, next_winding + 1] >= 0 && continue
                distance = distances[node, winding + 1] + 1
                distances[neighbor, next_winding + 1] = distance
                push!(queue, (neighbor, next_winding))
                neighbor == start && next_winding != 0 &&
                    (shortest = min(shortest, distance))
            end
        end
    end
    return shortest
end

if isdefined(TopoNoise, :Toric3DEncodedSlice)
    @testset "Honeycomb seed and encoded 3D slab algebra" begin
        lattice = TopoNoise.Toric3DEncodedSlice
        @test_throws ArgumentError lattice(0, 2, 2)
        @test_throws ArgumentError lattice(2, 1, 2)
        for (t, ly, lz) in ((1, 2, 2), (2, 2, 2), (2, 2, 3), (2, 3, 3))
            code = lattice(t, ly, lz)
            m = ly * lz
            n = 3(t + 1)m
            @test data_qubit_count(code) == n
            @test all(data_qubit_index(code, data_qubit_coordinate(code, q)...) == q
                      for q in 1:n)
            @test data_qubit_coordinates(code) == data_qubit_coordinates(code.lattice)
            hx = star_check_matrix(code)
            hz = plaquette_check_matrix(code)
            @test size(hx) == ((t + 2)m, n)
            @test size(hz) == (3t*m, n)
            @test all(==(2), vec(sum(hx; dims=1)))
            @test all(iszero, Int.(hx) * transpose(Int.(hz)) .% 2)
            rx = toric3d_encoded_test_rank(hx)
            rz = toric3d_encoded_test_rank(hz)
            @test rx == (t + 2)m - 1
            @test rz == (2t + 1)m - 1
            @test n - rx - rz == 2
            @test Matrix(star_check_matrix(code; sparse=true)) == hx
            @test Matrix(plaquette_check_matrix(code; sparse=true)) == hz
            lx = TopoNoise._check_matrix(TopoNoise.toric3d_encoded_logical_x(code), n; sparse=false)
            lzop = TopoNoise._check_matrix(TopoNoise.toric3d_encoded_logical_z(code), n; sparse=false)
            @test all(iszero, Int.(lx) * transpose(Int.(hz)) .% 2)
            @test all(iszero, Int.(lzop) * transpose(Int.(hx)) .% 2)
            @test Int.(lx) * transpose(Int.(lzop)) .% 2 == [1 0; 0 1]
            @test toric3d_encoded_test_rank(vcat(hx, lx)) == rx + 2
            @test toric3d_encoded_test_rank(vcat(hz, lzop)) == rz + 2
            @test vec(sum(lx; dims=2)) == [(t + 1)lz, (t + 1)ly]
            @test vec(sum(lzop; dims=2)) == [2ly, 2lz]
            @test toric3d_encoded_test_z_distance(code) == 2min(ly, lz)

            seed = TopoNoise.toric3d_encoded_seed_checks(code)
            sx = TopoNoise._check_matrix(seed.x, 3m; sparse=false)
            sz = TopoNoise._check_matrix(seed.z, 3m; sparse=false)
            @test length(seed.x) == 2m
            @test length(seed.z) == m
            @test all(length(support) == 3 for support in seed.x)
            @test all(length(support) == 6 for support in seed.z)
            @test all(iszero, Int.(sx) * transpose(Int.(sz)) .% 2)
            @test toric3d_encoded_test_rank(sx) == 2m - 1
            @test toric3d_encoded_test_rank(sz) == m - 1
            @test 3m - toric3d_encoded_test_rank(sx) - toric3d_encoded_test_rank(sz) == 2
        end
    end

    @testset "Encoded transfer prepares logical plus at every checkpoint" begin
        for dims in ((1, 2, 2), (2, 2, 3), (3, 3, 3))
            code = TopoNoise.Toric3DEncodedSlice(dims...)
            encoder = TopoNoise.toric3d_encoded_encoder(code)
            @test encoder.code === code.lattice
            @test TopoNoise.verify_toric3d_encoded_encoder(code)
            @test TopoNoise.verify_toric3d_encoded_encoder(code, encoder)
            short = TopoNoise.Toric3DSliceEncoder(code.lattice,
                encoder.operations[1:end-1], encoder.checkpoints)
            @test !TopoNoise.verify_toric3d_encoded_encoder(code, short)
            checkpoints = TopoNoise.Toric3DSliceEncoder(code.lattice,
                encoder.operations, encoder.checkpoints .+ 1)
            @test !TopoNoise.verify_toric3d_encoded_encoder(code, checkpoints)
            changed = copy(encoder.operations)
            changed[1] = TopoNoise.EncoderOperation(:CNOT, changed[2].qubits)
            corrupted = TopoNoise.Toric3DSliceEncoder(code.lattice, changed, encoder.checkpoints)
            @test !TopoNoise.verify_toric3d_encoded_encoder(code, corrupted)
        end
    end
end
