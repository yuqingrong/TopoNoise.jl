"""
Diagonal cubic-lattice preparation with slices `0:steps`, periodic transverse
coordinates, and an initially fixed all-zero slice. These boundary conditions
specify a unique state; they do not encode the old rough slab's logical qubit.
"""
struct Toric3DSliceLattice
    steps::Int
    Ly::Int
    Lz::Int
    function Toric3DSliceLattice(steps::Integer, Ly::Integer, Lz::Integer)
        1 <= steps <= typemax(Int) ||
            throw(ArgumentError("steps must be at least 1 and fit Int"))
        2 <= Ly <= typemax(Int) && 2 <= Lz <= typemax(Int) ||
            throw(ArgumentError("periodic dimensions must be at least 2 and fit Int"))
        3 * big(Ly) * Lz * (big(steps) + 1) <= typemax(Int) ||
            throw(ArgumentError("the lattice qubit count must fit Int"))
        return new(Int(steps), Int(Ly), Int(Lz))
    end
end

const _TORIC3D_SLICE_SHIFTS = ((1, 0), (0, 1), (0, 0))

data_qubit_count(code::Toric3DSliceLattice) =
    3 * code.Ly * code.Lz * (code.steps + 1)

"""Index an incoming oriented edge ending at `(t,u,v)`, with direction `1:3`."""
function data_qubit_index(code::Toric3DSliceLattice, direction::Integer,
                          t::Integer, u::Integer, v::Integer)
    1 <= direction <= 3 && 0 <= t <= code.steps &&
        0 <= u < code.Ly && 0 <= v < code.Lz ||
        throw(ArgumentError("edge direction, slice, or periodic coordinate is out of range"))
    return (3 * Int(t) + Int(direction) - 1) * code.Ly * code.Lz +
           Int(u) * code.Lz + Int(v) + 1
end

"""Return `(direction,t,u,v)` for one diagonal-lattice data qubit."""
function data_qubit_coordinate(code::Toric3DSliceLattice, index::Integer)
    1 <= index <= data_qubit_count(code) ||
        throw(ArgumentError("data-qubit index is out of range"))
    block, offset = divrem(Int(index) - 1, code.Ly * code.Lz)
    t, direction = divrem(block, 3)
    u, v = divrem(offset, code.Lz)
    return (direction + 1, t, u, v)
end

data_qubit_coordinates(code::Toric3DSliceLattice) =
    [data_qubit_coordinate(code, q) for q in 1:data_qubit_count(code)]

"""Full data-qubit indices in one slice, ordered direction, u, then v."""
function toric3d_slice_qubits(code::Toric3DSliceLattice, t::Integer)
    0 <= t <= code.steps || throw(ArgumentError("slice is out of range"))
    return [data_qubit_index(code, direction, t, u, v)
            for direction in 1:3 for u in 0:code.Ly-1 for v in 0:code.Lz-1]
end

"""The old-edge controls of pair parity Kij for a new vertex `(t,u,v)`."""
function _toric3d_slice_pair(code::Toric3DSliceLattice, t::Integer,
                            u::Integer, v::Integer, i::Integer, j::Integer)
    1 <= t <= code.steps && 0 <= u < code.Ly && 0 <= v < code.Lz &&
        1 <= i <= 3 && 1 <= j <= 3 && i != j ||
        throw(ArgumentError("new vertex or distinct pair directions are out of range"))
    si, sj = _TORIC3D_SLICE_SHIFTS[i], _TORIC3D_SLICE_SHIFTS[j]
    return (data_qubit_index(code, j, t - 1, mod(u - si[1], code.Ly), mod(v - si[2], code.Lz)),
            data_qubit_index(code, i, t - 1, mod(u - sj[1], code.Ly), mod(v - sj[2], code.Lz)))
end

"""Bulk six-edge X stars, followed by three-edge stars on the final frontier."""
function star_checks(code::Toric3DSliceLattice)
    checks = Vector{Vector{Int}}()
    for t in 1:code.steps, u in 0:code.Ly-1, v in 0:code.Lz-1
        support = [data_qubit_index(code, i, t, u, v) for i in 1:3]
        if t < code.steps
            for (i, shift) in enumerate(_TORIC3D_SLICE_SHIFTS)
                push!(support, data_qubit_index(code, i, t + 1,
                    mod(u + shift[1], code.Ly), mod(v + shift[2], code.Lz)))
            end
        end
        push!(checks, support)
    end
    return checks
end

"""Four-edge physical plaquettes, ordered time, u, v, then pair 12,13,23."""
function plaquette_checks(code::Toric3DSliceLattice)
    checks = Vector{Vector{Int}}()
    for t in 1:code.steps, u in 0:code.Ly-1, v in 0:code.Lz-1
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            controls = _toric3d_slice_pair(code, t, u, v, i, j)
            push!(checks, [data_qubit_index(code, i, t, u, v),
                          data_qubit_index(code, j, t, u, v), controls...])
        end
    end
    return checks
end

star_check_matrix(code::Toric3DSliceLattice; sparse::Bool=false) =
    _check_matrix(star_checks(code), data_qubit_count(code); sparse)
plaquette_check_matrix(code::Toric3DSliceLattice; sparse::Bool=false) =
    _check_matrix(plaquette_checks(code), data_qubit_count(code); sparse)

"""Individual Z checks fixing the initial product-state boundary."""
toric3d_slice_boundary_checks(code::Toric3DSliceLattice) =
    [[q] for q in toric3d_slice_qubits(code, 0)]

"""
Six-edge Z compatibility loops contained in slice `t`, ordered by `(u,v)`.
Each physical edge occurs in two checks, including distinct periodic edges at
transverse length two. Together these checks form a triangular matching graph.
"""
function toric3d_slice_checks(code::Toric3DSliceLattice, t::Integer)
    0 <= t <= code.steps || throw(ArgumentError("slice is out of range"))
    checks = Vector{Vector{Int}}()
    for u in 0:code.Ly-1, v in 0:code.Lz-1
        support = Int[]
        for (i, j) in ((1, 2), (1, 3), (2, 3))
            si, sj = _TORIC3D_SLICE_SHIFTS[i], _TORIC3D_SLICE_SHIFTS[j]
            push!(support,
                data_qubit_index(code, j, t, mod(u - si[1], code.Ly), mod(v - si[2], code.Lz)),
                data_qubit_index(code, i, t, mod(u - sj[1], code.Ly), mod(v - sj[2], code.Lz)))
        end
        push!(checks, support)
    end
    return checks
end

"""The local N-by-3N compatibility matrix, with columns direction,u,v."""
toric3d_slice_check_matrix(code::Toric3DSliceLattice; sparse::Bool=false) =
    _check_matrix(toric3d_slice_checks(code, 0), 3 * code.Ly * code.Lz; sparse)
