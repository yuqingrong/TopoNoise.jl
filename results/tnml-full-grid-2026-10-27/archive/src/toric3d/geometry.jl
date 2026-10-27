"""Relative cubic toric-code slab with rough x ends and periodic y/z."""
struct ToricCode3D
    Lx::Int
    Ly::Int
    Lz::Int
    function ToricCode3D(Lx::Integer, Ly::Integer, Lz::Integer)
        1 <= Lx <= typemax(Int) ||
            throw(ArgumentError("Lx must be at least 1 and fit Int"))
        2 <= Ly <= typemax(Int) ||
            throw(ArgumentError("Ly must be at least 2 and fit Int"))
        2 <= Lz <= typemax(Int) ||
            throw(ArgumentError("Lz must be at least 2 and fit Int"))
        return new(Int(Lx), Int(Ly), Int(Lz))
    end
end

data_qubit_count(code::ToricCode3D) = (3 * code.Lx + 1) * code.Ly * code.Lz

"""Return the index of an oriented edge at its lower/starting coordinate."""
function data_qubit_index(code::ToricCode3D, direction::Symbol,
                          x::Integer, y::Integer, z::Integer)
    0 <= y < code.Ly && 0 <= z < code.Lz ||
        throw(ArgumentError("periodic coordinates are out of range"))
    plane = code.Ly * code.Lz
    offset = Int(y) * code.Lz + Int(z) + 1
    if direction === :x && 0 <= x <= code.Lx
        return Int(x) * plane + offset
    elseif direction === :y && 1 <= x <= code.Lx
        return (code.Lx + Int(x)) * plane + offset
    elseif direction === :z && 1 <= x <= code.Lx
        return (2 * code.Lx + Int(x)) * plane + offset
    end
    throw(ArgumentError("edge orientation or x coordinate is out of range"))
end

"""Return `(direction, x, y, z)` for an oriented edge index."""
function data_qubit_coordinate(code::ToricCode3D, index::Integer)
    1 <= index <= data_qubit_count(code) ||
        throw(ArgumentError("data-qubit index is out of range"))
    plane = code.Ly * code.Lz
    block, offset = divrem(Int(index) - 1, plane)
    y, z = divrem(offset, code.Lz)
    if block <= code.Lx
        return (:x, block, y, z)
    elseif block <= 2 * code.Lx
        return (:y, block - code.Lx, y, z)
    end
    return (:z, block - 2 * code.Lx, y, z)
end

data_qubit_coordinates(code::ToricCode3D) =
    [data_qubit_coordinate(code, q) for q in 1:data_qubit_count(code)]

"""X stars at all retained vertices, ordered x then y then z."""
function star_checks(code::ToricCode3D)
    q(d, x, y, z) = data_qubit_index(code, d, x, mod(y, code.Ly), mod(z, code.Lz))
    return [Int[
        q(:x, x - 1, y, z), q(:x, x, y, z),
        q(:y, x, y, z), q(:y, x, y - 1, z),
        q(:z, x, y, z), q(:z, x, y, z - 1),
    ] for x in 1:code.Lx for y in 0:code.Ly-1 for z in 0:code.Lz-1]
end

"""Z plaquettes, including triangular cells on both rough x ends."""
function plaquette_checks(code::ToricCode3D)
    q(d, x, y, z) = data_qubit_index(code, d, x, mod(y, code.Ly), mod(z, code.Lz))
    checks = Vector{Vector{Int}}()
    for x in 0:code.Lx, y in 0:code.Ly-1, z in 0:code.Lz-1
        support = Int[q(:x, x, y, z), q(:x, x, y + 1, z)]
        x > 0 && push!(support, q(:y, x, y, z))
        x < code.Lx && push!(support, q(:y, x + 1, y, z))
        push!(checks, support)
    end
    for x in 0:code.Lx, y in 0:code.Ly-1, z in 0:code.Lz-1
        support = Int[q(:x, x, y, z), q(:x, x, y, z + 1)]
        x > 0 && push!(support, q(:z, x, y, z))
        x < code.Lx && push!(support, q(:z, x + 1, y, z))
        push!(checks, support)
    end
    for x in 1:code.Lx, y in 0:code.Ly-1, z in 0:code.Lz-1
        push!(checks, Int[
            q(:y, x, y, z), q(:y, x, y, z + 1),
            q(:z, x, y, z), q(:z, x, y + 1, z),
        ])
    end
    return checks
end

star_check_matrix(code::ToricCode3D; sparse::Bool=false) =
    _check_matrix(star_checks(code), data_qubit_count(code); sparse)
plaquette_check_matrix(code::ToricCode3D; sparse::Bool=false) =
    _check_matrix(plaquette_checks(code), data_qubit_count(code); sparse)

"""A minimum transverse X membrane across one x-edge plane."""
logical_x_support(code::ToricCode3D) =
    [data_qubit_index(code, :x, 0, y, z)
     for y in 0:code.Ly-1 for z in 0:code.Lz-1]

"""A minimum longitudinal Z string from one rough end to the other."""
logical_z_support(code::ToricCode3D) =
    [data_qubit_index(code, :x, x, 0, 0) for x in 0:code.Lx]

"""The prepared prefix code and its local-to-full edge map after layer `j`."""
function toric3d_prefix(code::ToricCode3D, j::Integer)
    1 <= j <= code.Lx || throw(ArgumentError("prefix layer is out of range"))
    prefix = ToricCode3D(j, code.Ly, code.Lz)
    qubits = [data_qubit_index(code, data_qubit_coordinate(prefix, q)...)
              for q in 1:data_qubit_count(prefix)]
    return (code=prefix, qubits=qubits)
end
