"""A CSS stabilizer check with its conventional name, Pauli type, and support."""
struct StabilizerCheck
    name::Symbol
    pauli::Symbol
    support::Vector{Int}
end

"""
    RotatedPlanarCode(distance; boundary_orientation=:x_ns)

Geometry for an odd-distance rotated planar surface code.  `:x_ns` has
Z-type A_s boundary checks on the north and south edges; `:x_ew` is its
90-degree clockwise rotation.
"""
struct RotatedPlanarCode
    distance::Int
    boundary_orientation::Symbol
    a_s::Vector{Vector{Int}}
    b_p::Vector{Vector{Int}}
end

function RotatedPlanarCode(distance::Integer; boundary_orientation::Symbol=:x_ns)
    distance < 3 && throw(ArgumentError("distance must be an odd integer at least 3"))
    iseven(distance) && throw(ArgumentError("distance must be an odd integer at least 3"))
    distance > typemax(Int) && throw(ArgumentError("distance is unsupported on this platform"))
    boundary_orientation in (:x_ns, :x_ew) ||
        throw(ArgumentError("boundary_orientation must be :x_ns or :x_ew"))

    d = Int(distance)
    a_s, b_p = _x_ns_checks(d)
    if boundary_orientation === :x_ew
        a_s = [_rotate_clockwise(support, d) for support in a_s]
        b_p = [_rotate_clockwise(support, d) for support in b_p]
    end
    return RotatedPlanarCode(d, boundary_orientation, a_s, b_p)
end

distance(code::RotatedPlanarCode) = code.distance
boundary_orientation(code::RotatedPlanarCode) = code.boundary_orientation
data_qubit_count(code::RotatedPlanarCode) = code.distance^2

function data_qubit_index(code::RotatedPlanarCode, row::Integer, column::Integer)
    d = distance(code)
    1 <= row <= d && 1 <= column <= d ||
        throw(ArgumentError("data-qubit coordinates must lie in 1:$d by 1:$d"))
    return (Int(row) - 1) * d + Int(column)
end

function data_qubit_coordinate(code::RotatedPlanarCode, index::Integer)
    1 <= index <= data_qubit_count(code) ||
        throw(ArgumentError("data-qubit index must lie in 1:$(data_qubit_count(code))"))
    d = distance(code)
    return (fld(Int(index) - 1, d) + 1, mod(Int(index) - 1, d) + 1)
end

data_qubit_coordinates(code::RotatedPlanarCode) =
    [(row, column) for row in 1:distance(code) for column in 1:distance(code)]

a_s_checks(code::RotatedPlanarCode) = copy.(code.a_s)
b_p_checks(code::RotatedPlanarCode) = copy.(code.b_p)

function stabilizers(code::RotatedPlanarCode)
    return vcat(
        [StabilizerCheck(:A_s, :Z, copy(support)) for support in code.a_s],
        [StabilizerCheck(:B_p, :X, copy(support)) for support in code.b_p],
    )
end

a_s_check_matrix(code::RotatedPlanarCode; sparse::Bool=false) =
    _check_matrix(code.a_s, data_qubit_count(code); sparse)
b_p_check_matrix(code::RotatedPlanarCode; sparse::Bool=false) =
    _check_matrix(code.b_p, data_qubit_count(code); sparse)
stabilizer_check_matrix(code::RotatedPlanarCode; sparse::Bool=false) =
    _check_matrix(vcat(code.a_s, code.b_p), data_qubit_count(code); sparse)

function logical_x_support(code::RotatedPlanarCode)
    d = distance(code)
    if boundary_orientation(code) === :x_ns
        return [data_qubit_index(code, 1, column) for column in 1:d]
    end
    return [data_qubit_index(code, row, d) for row in 1:d]
end

function logical_z_support(code::RotatedPlanarCode)
    d = distance(code)
    if boundary_orientation(code) === :x_ns
        return [data_qubit_index(code, row, 1) for row in 1:d]
    end
    return [data_qubit_index(code, 1, column) for column in 1:d]
end

function _x_ns_checks(d::Int)
    a_s = Vector{Vector{Int}}()
    b_p = Vector{Vector{Int}}()
    q(row, column) = (row - 1) * d + column

    for row in 1:d-1, column in 1:d-1
        support = [q(row, column), q(row, column + 1),
                   q(row + 1, column), q(row + 1, column + 1)]
        if iseven(row + column)
            push!(a_s, support)
        else
            push!(b_p, support)
        end
    end
    for column in 1:d-1
        isodd(1 + column) && push!(a_s, [q(1, column), q(1, column + 1)])
    end
    for column in 1:d-1
        isodd(d - 1 + column) && push!(a_s, [q(d, column), q(d, column + 1)])
    end
    for row in 1:d-1
        iseven(row + 1) && push!(b_p, [q(row, 1), q(row + 1, 1)])
    end
    for row in 1:d-1
        iseven(row + d - 1) && push!(b_p, [q(row, d), q(row + 1, d)])
    end
    return a_s, b_p
end

function _rotate_clockwise(support::Vector{Int}, d::Int)
    rotated = Int[]
    for index in support
        row = fld(index - 1, d) + 1
        column = mod(index - 1, d) + 1
        push!(rotated, (column - 1) * d + (d + 1 - row))
    end
    sort!(rotated)
    return rotated
end

function _check_matrix(checks::Vector{Vector{Int}}, qubit_count::Int; sparse::Bool)
    matrix = falses(length(checks), qubit_count)
    for (row, support) in enumerate(checks), qubit in support
        matrix[row, qubit] = true
    end
    return sparse ? SparseArrays.sparse(matrix) : matrix
end
