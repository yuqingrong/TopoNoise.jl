using Test
using TopoNoise

function toric3d_test_rank(matrix)
    work = BitMatrix(matrix)
    rows, columns = size(work)
    pivot = 0
    for column in 1:columns
        candidate = findfirst(row -> work[row, column], pivot+1:rows)
        isnothing(candidate) && continue
        source = pivot + candidate
        pivot += 1
        work[pivot, :], work[source, :] = copy(work[source, :]), copy(work[pivot, :])
        for row in pivot+1:rows
            work[row, column] && (work[row, :] .⊻= work[pivot, :])
        end
    end
    return pivot
end

@testset "3D relative toric slab geometry" begin
    @test_throws ArgumentError ToricCode3D(0, 2, 2)
    @test_throws ArgumentError ToricCode3D(1, 1, 2)
    @test_throws ArgumentError ToricCode3D(1, 2, 1)
    for (lx, ly, lz) in ((1, 2, 2), (2, 2, 2), (2, 3, 2))
        code = ToricCode3D(lx, ly, lz)
        n = (3lx + 1) * ly * lz
        @test data_qubit_count(code) == n
        @test length(star_checks(code)) == lx * ly * lz
        @test length(plaquette_checks(code)) == (3lx + 2) * ly * lz
        @test all(length(check) == 6 for check in star_checks(code))
        @test count(check -> length(check) == 3, plaquette_checks(code)) ==
              4 * ly * lz
        @test all(length(check) in (3, 4) for check in plaquette_checks(code))
        hx = star_check_matrix(code)
        hz = plaquette_check_matrix(code)
        @test size(hx) == (lx * ly * lz, n)
        @test size(hz) == ((3lx + 2) * ly * lz, n)
        @test all(iszero, (Int.(hx) * transpose(Int.(hz))) .% 2)
        @test toric3d_test_rank(hx) == lx * ly * lz
        @test toric3d_test_rank(hz) == n - lx * ly * lz - 1
        logical_x = logical_x_support(code)
        logical_z = logical_z_support(code)
        @test length(logical_x) == ly * lz
        @test length(logical_z) == lx + 1
        @test isodd(length(intersect(logical_x, logical_z)))
        @test all(iszero, (Int.(hz) * Int.([q in logical_x for q in 1:n])) .% 2)
        @test all(iszero, (Int.(hx) * Int.([q in logical_z for q in 1:n])) .% 2)
        @test length(unique(data_qubit_coordinate(code, q) for q in 1:n)) == n
        @test all(data_qubit_index(code, data_qubit_coordinate(code, q)...) == q
                  for q in 1:n)
    end
    code = ToricCode3D(2, 2, 2)
    @test data_qubit_count(code) == 28
    @test length(plaquette_checks(code)) == 32
    for j in 1:2
        prefix = toric3d_prefix(code, j)
        @test prefix.code.Lx == j
        @test length(prefix.qubits) == data_qubit_count(prefix.code)
        @test all(data_qubit_coordinate(prefix.code, q) ==
                  data_qubit_coordinate(code, prefix.qubits[q])
                  for q in eachindex(prefix.qubits))
    end
end
