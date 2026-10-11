using Test
using TopoNoise

function toric3d_slice_test_rank(matrix)
    work = BitMatrix(matrix)
    rows, columns = size(work)
    pivot = 0
    for column in 1:columns
        source = findfirst(row -> work[row, column], pivot+1:rows)
        isnothing(source) && continue
        source += pivot
        pivot += 1
        work[pivot, :], work[source, :] = copy(work[source, :]), copy(work[pivot, :])
        for row in pivot+1:rows
            work[row, column] && (work[row, :] .⊻= work[pivot, :])
        end
    end
    return pivot
end

@testset "Diagonal 3D slice geometry and physical boundary state" begin
    @test isdefined(TopoNoise, :Toric3DSliceLattice)
    if isdefined(TopoNoise, :Toric3DSliceLattice)
        lattice = TopoNoise.Toric3DSliceLattice
        @test_throws ArgumentError lattice(0, 2, 2)
        @test_throws ArgumentError lattice(2, 1, 2)
        @test_throws ArgumentError lattice(2, 2, 1)
        for (steps, ly, lz) in ((1, 2, 2), (2, 2, 2), (3, 2, 3), (2, 3, 3))
            code = lattice(steps, ly, lz)
            cells = ly * lz
            n = 3 * cells * (steps + 1)
            @test data_qubit_count(code) == n
            @test all(data_qubit_index(code, data_qubit_coordinate(code, q)...) == q
                      for q in 1:n)
            @test length(unique(data_qubit_coordinate(code, q) for q in 1:n)) == n
            @test vcat([TopoNoise.toric3d_slice_qubits(code, t) for t in 0:steps]...) == 1:n
            @test_throws ArgumentError data_qubit_index(code, 0, 0, 0, 0)
            @test_throws ArgumentError data_qubit_index(code, 4, 0, 0, 0)
            @test_throws ArgumentError data_qubit_index(code, 1, -1, 0, 0)
            @test_throws ArgumentError data_qubit_index(code, 1, steps + 1, 0, 0)
            @test_throws ArgumentError data_qubit_index(code, 1, 0, ly, 0)
            @test_throws ArgumentError data_qubit_index(code, 1, 0, 0, -1)
            @test_throws ArgumentError data_qubit_coordinate(code, 0)
            @test_throws ArgumentError data_qubit_coordinate(code, n + 1)
            @test_throws ArgumentError TopoNoise.toric3d_slice_qubits(code, steps + 1)
            @test_throws ArgumentError TopoNoise.toric3d_slice_checks(code, -1)
            @test_throws ArgumentError TopoNoise._toric3d_slice_pair(code, 0, 0, 0, 1, 2)
            @test_throws ArgumentError TopoNoise._toric3d_slice_pair(code, 1, 0, 0, 2, 2)

            stars = star_checks(code)
            plaquettes = plaquette_checks(code)
            boundary = TopoNoise.toric3d_slice_boundary_checks(code)
            @test length(stars) == steps * cells
            @test count(check -> length(check) == 3, stars) == cells
            @test count(check -> length(check) == 6, stars) == (steps - 1) * cells
            @test length(plaquettes) == 3 * steps * cells
            @test all(length(check) == length(unique(check)) == 4 for check in plaquettes)
            @test boundary == [[q] for q in 1:3cells]
            hx = star_check_matrix(code)
            hp = plaquette_check_matrix(code)
            hz = TopoNoise._check_matrix(vcat(plaquettes, boundary), n; sparse=false)
            @test all(iszero, (Int.(hx) * transpose(Int.(hz))) .% 2)
            @test toric3d_slice_test_rank(hx) == steps * cells
            @test toric3d_slice_test_rank(hp) == (2steps + 1) * cells - 1
            @test toric3d_slice_test_rank(hz) == (2steps + 3) * cells
            @test toric3d_slice_test_rank(hx) + toric3d_slice_test_rank(hz) == n
            @test Matrix(star_check_matrix(code; sparse=true)) == hx
            @test Matrix(plaquette_check_matrix(code; sparse=true)) == hp

            q = TopoNoise.toric3d_slice_check_matrix(code)
            @test size(q) == (cells, 3cells)
            @test all(==(2), vec(sum(q; dims=1)))
            @test all(==(6), vec(sum(q; dims=2)))
            @test toric3d_slice_test_rank(q) == cells - 1
            for t in 0:steps
                checks = TopoNoise.toric3d_slice_checks(code, t)
                @test all(length(check) == length(unique(check)) == 6 for check in checks)
                @test all(all(data_qubit_coordinate(code, id)[2] == t for id in check)
                          for check in checks)
                hq = TopoNoise._check_matrix(checks, n; sparse=false)
                @test toric3d_slice_test_rank(vcat(hz, hq)) == toric3d_slice_test_rank(hz)
            end
        end

        small = lattice(2, 2, 2)
        @test data_qubit_count(small) == 36
        @test TopoNoise._toric3d_slice_pair(small, 1, 0, 0, 1, 2) == (7, 2)
        @test first(TopoNoise.toric3d_slice_checks(small, 0)) == [7, 2, 11, 1, 10, 5]
        @test first(star_checks(small)) == [13, 17, 21, 27, 30, 33]
        @test star_checks(small)[5] == [25, 29, 33]
        @test plaquette_checks(small)[1:3] == [[13, 17, 7, 2], [13, 21, 11, 1], [17, 21, 10, 5]]
        q = TopoNoise.toric3d_slice_check_matrix(small)
        @test size(q, 2) == 12
        @test length(unique(Tuple(q[:, column]) for column in axes(q, 2))) == 6
    end
end
