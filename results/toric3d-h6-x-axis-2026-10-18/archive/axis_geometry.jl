# Included inside AxisH6. These are new methods for a new geometry type; the
# historical diagonal and rough-boundary types and schedules remain unchanged.
import TopoNoise: star_checks, plaquette_checks

"""
Periodic yz square-toric seed followed by T physical-x growth transfers.
Q0 contains 2LyLz seed edges; each subsequent Qr contains x, y, z edge blocks.
The two transverse logical qubits persist for every prefix, including T=0.
"""
struct AxisCode
    T::Int
    Ly::Int
    Lz::Int
    function AxisCode(T::Integer, Ly::Integer, Lz::Integer)
        T >= 0 || throw(ArgumentError("T must be nonnegative"))
        Ly >= 2 || throw(ArgumentError("Ly must be at least two"))
        Lz >= 2 || throw(ArgumentError("Lz must be at least two"))
        new(Int(T), Int(Ly), Int(Lz))
    end
end

"""Literal Clifford operations and a checkpoint after each complete transfer."""
struct AxisEncoder
    code::AxisCode
    operations::Vector{TN.EncoderOperation}
    checkpoints::Vector{Int}
end

nqubits(code::AxisCode) = (3code.T + 2) * code.Ly * code.Lz

"""
Index a physical edge with periodic u,v coordinates. An x edge at r joins
vertices (r-1,u,v) and (r,u,v); y/z edges at r lie in that yz plane.
Within each direction block u is outermost and v is innermost.
"""
function edge(code::AxisCode, direction::Symbol, r::Integer, u::Integer, v::Integer)
    direction in (:x, :y, :z) || throw(ArgumentError("edge direction must be :x, :y, or :z"))
    minimum_r = direction === :x ? 1 : 0
    minimum_r <= r <= code.T || throw(ArgumentError("edge plane outside the code"))
    N = code.Ly * code.Lz
    cell = mod(u, code.Ly) * code.Lz + mod(v, code.Lz) + 1
    if r == 0
        return (direction === :y ? 0 : N) + cell
    end
    base = 2N + 3N * (r-1)
    offset = direction === :x ? 0 : direction === :y ? N : 2N
    return base + offset + cell
end

"""All y edges followed by all z edges of plane r, each in u,v order."""
function plane_edges(code::AxisCode, r::Integer)
    0 <= r <= code.T || throw(ArgumentError("plane outside the code"))
    return [edge(code, direction, r, u, v)
        for direction in (:y, :z) for u in 0:code.Ly-1 for v in 0:code.Lz-1]
end

"""Contiguous Q0 seed or Qr fresh x/y/z block, in physical-qubit index order."""
function slice_qubits(code::AxisCode, r::Integer)
    0 <= r <= code.T || throw(ArgumentError("slice outside the code"))
    N = code.Ly * code.Lz
    r == 0 && return collect(1:2N)
    return collect(2N + 3N*(r-1) + 1 : 2N + 3N*r)
end

function _axis_incident_plane_edges(code::AxisCode, r, u, v)
    return [edge(code, :y, r, u, v), edge(code, :y, r, u-1, v),
        edge(code, :z, r, u, v), edge(code, :z, r, u, v-1)]
end

"""
Physical X stars in plane-major r,u,v order. The seed has four-edge stars;
grown slabs have five-edge boundary stars and six-edge interior stars.
"""
function star_checks(code::AxisCode)
    checks = Vector{Int}[]
    for r in 0:code.T, u in 0:code.Ly-1, v in 0:code.Lz-1
        support = _axis_incident_plane_edges(code, r, u, v)
        r > 0 && push!(support, edge(code, :x, r, u, v))
        r < code.T && push!(support, edge(code, :x, r+1, u, v))
        push!(checks, support)
    end
    return checks
end

"""
Physical square Z plaquettes: all yz planes first, then all xy transfers,
then all xz transfers; within each family use r,u,v order. Period-two
parallel edges retain separate qubit indices, even when endpoints coincide.
"""
function plaquette_checks(code::AxisCode)
    checks = Vector{Int}[]
    for r in 0:code.T, u in 0:code.Ly-1, v in 0:code.Lz-1
        push!(checks, [edge(code, :y, r, u, v), edge(code, :y, r, u, v+1),
            edge(code, :z, r, u, v), edge(code, :z, r, u+1, v)])
    end
    for r in 1:code.T, u in 0:code.Ly-1, v in 0:code.Lz-1
        push!(checks, [edge(code, :y, r-1, u, v), edge(code, :y, r, u, v),
            edge(code, :x, r, u, v), edge(code, :x, r, u+1, v)])
    end
    for r in 1:code.T, u in 0:code.Ly-1, v in 0:code.Lz-1
        push!(checks, [edge(code, :z, r-1, u, v), edge(code, :z, r, u, v),
            edge(code, :x, r, u, v), edge(code, :x, r, u, v+1)])
    end
    return checks
end

"""Transverse dual X membranes, paired with the y and z seed Z cycles."""
function logical_x(code::AxisCode)
    return [[edge(code, :y, r, 0, v) for r in 0:code.T for v in 0:code.Lz-1],
        [edge(code, :z, r, u, 0) for r in 0:code.T for u in 0:code.Ly-1]]
end

"""Two Z cycles on the seed plane winding once around y and z."""
function logical_z(code::AxisCode)
    return [[edge(code, :y, 0, u, 0) for u in 0:code.Ly-1],
        [edge(code, :z, 0, 0, v) for v in 0:code.Lz-1]]
end

"""
The physical-x H+6 schedule. First copy every old y/z edge to its fresh
counterpart, in plane_edges order. Only after all copies, visit vertices in
u,v order: H on the fresh x link, then CNOTs to +y,-y,+z,-z fresh edges.
Every transfer has N H gates, 6N CNOTs and 13N independent gate outputs.
No H or CNOT target belongs to an earlier slice.
"""
function axis_encoder(code::AxisCode)
    operations = TN.EncoderOperation[]
    checkpoints = Int[]
    for r in 1:code.T
        for (old, fresh) in zip(plane_edges(code, r-1), plane_edges(code, r))
            push!(operations, TN.EncoderOperation(:CNOT, [old, fresh]))
        end
        for u in 0:code.Ly-1, v in 0:code.Lz-1
            link = edge(code, :x, r, u, v)
            push!(operations, TN.EncoderOperation(:H, [link]))
            for target in _axis_incident_plane_edges(code, r, u, v)
                push!(operations, TN.EncoderOperation(:CNOT, [link, target]))
            end
        end
        push!(checkpoints, length(operations))
    end
    return AxisEncoder(code, operations, checkpoints)
end

"""
Verify every prefix's physical stabilizers and the prepared logical signs
by signed Clifford pullback. basis=:zero means a square-toric |00> seed;
basis=:plus means |++>. Every nonseed qubit starts in |0>. Also validate
gate arity/counts, prefix bounds and that every target is in the fresh slice.
This checks the stored circuit algebra, rather than trusting its generator.
"""
function verify_encoder(code::AxisCode, encoder::AxisEncoder=axis_encoder(code);
        basis::Symbol=:zero)
    basis in (:zero, :plus) || throw(ArgumentError("basis must be :zero or :plus"))
    actual = encoder.code
    (actual.T, actual.Ly, actual.Lz) == (code.T, code.Ly, code.Lz) || return false
    N = code.Ly * code.Lz
    encoder.checkpoints == [7N*r for r in 1:code.T] || return false
    length(encoder.operations) == 7N*code.T || return false
    n = nqubits(code)
    for operation in encoder.operations
        arity = operation.gate === :H ? 1 : operation.gate === :CNOT ? 2 : 0
        arity > 0 && length(operation.qubits) == arity || return false
        all(q -> 1 <= q <= n, operation.qubits) || return false
        length(unique(operation.qubits)) == arity || return false
    end
    for r in 1:code.T
        operations = @view encoder.operations[7N*(r-1)+1:7N*r]
        count(op -> op.gate === :H, operations) == N || return false
        count(op -> op.gate === :CNOT, operations) == 6N || return false
        firstfresh, lastfresh = first(slice_qubits(code,r)), last(slice_qubits(code,r))
        all(firstfresh <= last(op.qubits) <= lastfresh for op in operations) || return false
        all(all(q <= lastfresh for q in op.qubits) for op in operations) || return false
    end
    seed = AxisCode(0, code.Ly, code.Lz)
    seed_x_checks, seed_z_checks = star_checks(seed), plaquette_checks(seed)
    if basis === :zero
        append!(seed_z_checks, logical_z(seed))
    else
        append!(seed_x_checks, logical_x(seed))
    end
    seed_x = TN._check_matrix(seed_x_checks, 2N; sparse=false)
    seed_z = TN._check_matrix(seed_z_checks, 2N; sparse=false)
    seed_x_basis = TN._toric3d_encoded_row_basis(seed_x)
    seed_z_basis = TN._toric3d_encoded_row_basis(seed_z)
    logical_pauli = basis === :zero ? :Z : :X
    for r in 0:code.T
        prefix = AxisCode(r, code.Ly, code.Lz)
        operations = @view encoder.operations[1:7N*r]
        logical_checks = basis === :zero ? logical_z(prefix) : logical_x(prefix)
        for (pauli, checks) in ((:X, star_checks(prefix)), (:Z, plaquette_checks(prefix)),
                               (logical_pauli, logical_checks))
            all(TN._toric3d_encoded_pulled_to_seed(n, operations, pauli, support,
                2N, seed_x, seed_x_basis, seed_z_basis) for support in checks) || return false
        end
    end
    return true
end
