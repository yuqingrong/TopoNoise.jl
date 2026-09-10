using Test
using SparseArrays

function gf2_rank(matrix)
    rows = [BitVector(matrix[row, :]) for row in axes(matrix, 1)]
    pivot_row = 1
    for column in axes(matrix, 2)
        found = findfirst(row -> row[column], @view rows[pivot_row:end])
        isnothing(found) && continue
        pivot = pivot_row + found - 1
        rows[pivot], rows[pivot_row] = rows[pivot_row], rows[pivot]
        for row in eachindex(rows)
            row != pivot_row && rows[row][column] && (rows[row] .⊻= rows[pivot_row])
        end
        pivot_row += 1
        pivot_row > length(rows) && return length(rows)
    end
    return pivot_row - 1
end

@testset "Rotated planar code geometry" begin
    @testset "constructor validates odd distances and metadata" begin
        @test_throws ArgumentError RotatedPlanarCode(2)
        @test_throws ArgumentError RotatedPlanarCode(4)
        @test_throws ArgumentError RotatedPlanarCode(1)
        @test_throws ArgumentError RotatedPlanarCode(3; boundary_orientation=:diagonal)

        code = RotatedPlanarCode(3)
        @test distance(code) == 3
        @test boundary_orientation(code) == :x_ns
        @test data_qubit_count(code) == 9
    end

    @testset "d=3 :x_ns has hand-checked supports and check types" begin
        code = RotatedPlanarCode(3)

        # These supports are derived directly from the 3-by-3 fixture, not
        # from geometry helpers.
        @test a_s_checks(code) == [[2, 3, 5, 6], [4, 5, 7, 8], [6, 9], [1, 4]]
        @test b_p_checks(code) == [[5, 6, 8, 9], [1, 2, 4, 5], [2, 3], [7, 8]]
        @test [(check.name, check.pauli, check.support) for check in stabilizers(code)] == [
            (:A_s, :Z, [2, 3, 5, 6]), (:A_s, :Z, [4, 5, 7, 8]),
            (:A_s, :Z, [6, 9]), (:A_s, :Z, [1, 4]),
            (:B_p, :X, [5, 6, 8, 9]), (:B_p, :X, [1, 2, 4, 5]),
            (:B_p, :X, [2, 3]), (:B_p, :X, [7, 8]),
        ]
        @test filter(support -> length(support) == 2, b_p_checks(code)) == [[2, 3], [7, 8]]
        @test filter(support -> length(support) == 2, a_s_checks(code)) == [[6, 9], [1, 4]]
        @test sort(length.(a_s_checks(code))) == [2, 2, 4, 4]
        @test sort(length.(b_p_checks(code))) == [2, 2, 4, 4]
    end

    @testset "d=3 :x_ew has orthogonal hand-checked supports" begin
        code = RotatedPlanarCode(3; boundary_orientation=:x_ew)

        @test a_s_checks(code) == [[1, 2, 4, 5], [5, 6, 8, 9], [2, 3], [7, 8]]
        @test b_p_checks(code) == [[2, 3, 5, 6], [4, 5, 7, 8], [1, 4], [6, 9]]
        @test [(check.name, check.pauli) for check in stabilizers(code)] ==
            vcat(fill((:A_s, :Z), 4), fill((:B_p, :X), 4))
        @test filter(support -> length(support) == 2, b_p_checks(code)) == [[1, 4], [6, 9]]
        @test filter(support -> length(support) == 2, a_s_checks(code)) == [[2, 3], [7, 8]]
        @test sort(length.(a_s_checks(code))) == [2, 2, 4, 4]
        @test sort(length.(b_p_checks(code))) == [2, 2, 4, 4]
    end

    @testset "check counts and matrices scale with distance" begin
        for d in (3, 5, 7)
            code = RotatedPlanarCode(d)
            expected_per_type = (d^2 - 1) ÷ 2
            expected_boundary = d - 1
            expected_bulk = (d - 1)^2 ÷ 2

            @test length(a_s_checks(code)) == expected_per_type
            @test length(b_p_checks(code)) == expected_per_type
            @test length(stabilizers(code)) == d^2 - 1
            @test count(==(2), length.(a_s_checks(code))) == expected_boundary
            @test count(==(2), length.(b_p_checks(code))) == expected_boundary
            @test count(==(4), length.(a_s_checks(code))) == expected_bulk
            @test count(==(4), length.(b_p_checks(code))) == expected_bulk
            @test size(a_s_check_matrix(code)) == (expected_per_type, d^2)
            @test size(b_p_check_matrix(code)) == (expected_per_type, d^2)
            @test size(stabilizer_check_matrix(code)) == (d^2 - 1, d^2)
            @test eltype(a_s_check_matrix(code)) == Bool
            @test issparse(a_s_check_matrix(code; sparse=true))
        end
    end

    @testset "CSS checks commute and have full GF(2) rank" begin
        for d in (3, 5, 7), orientation in (:x_ns, :x_ew)
            code = RotatedPlanarCode(d; boundary_orientation=orientation)
            for a_support in a_s_checks(code), b_support in b_p_checks(code)
                @test iseven(length(intersect(a_support, b_support)))
            end
            @test gf2_rank(stabilizer_check_matrix(code)) == d^2 - 1
        end
    end

    @testset "canonical logical strings have the required commutation" begin
        logical_x_fixtures = Dict(:x_ns => [3, 6, 9], :x_ew => [1, 2, 3])
        logical_z_fixtures = Dict(:x_ns => [1, 2, 3], :x_ew => [1, 4, 7])
        for d in (3, 5, 7), orientation in (:x_ns, :x_ew)
            code = RotatedPlanarCode(d; boundary_orientation=orientation)
            logical_x = logical_x_support(code)
            logical_z = logical_z_support(code)

            @test length(logical_x) == d
            @test length(logical_z) == d
            if d == 3
                @test logical_x == logical_x_fixtures[orientation]
                @test logical_z == logical_z_fixtures[orientation]
            end
            @test isodd(length(intersect(logical_x, logical_z)))
            for a_support in a_s_checks(code)
                @test iseven(length(intersect(logical_x, a_support)))
            end
            for b_support in b_p_checks(code)
                @test iseven(length(intersect(logical_z, b_support)))
            end
        end
    end

    @testset "data-qubit indexing is row-major and invertible" begin
        code = RotatedPlanarCode(5)
        @test data_qubit_index(code, 1, 1) == 1
        @test data_qubit_index(code, 3, 4) == 14
        @test data_qubit_index(code, 5, 5) == 25
        @test data_qubit_coordinates(code) == [(row, col) for row in 1:5 for col in 1:5]
        for index in 1:data_qubit_count(code)
            row, col = data_qubit_coordinate(code, index)
            @test data_qubit_index(code, row, col) == index
        end
        @test_throws ArgumentError data_qubit_index(code, 0, 1)
        @test_throws ArgumentError data_qubit_coordinate(code, 26)
    end
end
