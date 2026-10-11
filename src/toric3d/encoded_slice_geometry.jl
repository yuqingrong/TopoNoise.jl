"""
Diagonal 3D transfer slab with a honeycomb toric-code seed boundary. It retains
two transverse logical qubits: `n=3(T+1)LyLz`, `dZ=2min(Ly,Lz)`, and
`dX=(T+1)min(Ly,Lz)`. The seed is prepared in logical `|++>` so a residual Z
logical operator is an actual state-preparation error.
"""
struct Toric3DEncodedSlice
    lattice::Toric3DSliceLattice
end

Toric3DEncodedSlice(steps::Integer, Ly::Integer, Lz::Integer) =
    Toric3DEncodedSlice(Toric3DSliceLattice(steps, Ly, Lz))

data_qubit_count(code::Toric3DEncodedSlice) = data_qubit_count(code.lattice)
data_qubit_index(code::Toric3DEncodedSlice, coordinates...) =
    data_qubit_index(code.lattice, coordinates...)
data_qubit_coordinate(code::Toric3DEncodedSlice, index::Integer) =
    data_qubit_coordinate(code.lattice, index)
data_qubit_coordinates(code::Toric3DEncodedSlice) = data_qubit_coordinates(code.lattice)

"""Honeycomb seed stars (black, then white) and six-edge Z hexagons."""
function toric3d_encoded_seed_checks(code::Toric3DEncodedSlice)
    lattice = code.lattice
    black, white = Vector{Int}[], Vector{Int}[]
    for u in 0:lattice.Ly-1, v in 0:lattice.Lz-1
        push!(black, [data_qubit_index(code, d, 0, u, v) for d in 1:3])
        push!(white, [data_qubit_index(code, d, 0,
            mod(u + shift[1], lattice.Ly), mod(v + shift[2], lattice.Lz))
            for (d, shift) in enumerate(_TORIC3D_SLICE_SHIFTS)])
    end
    return (x=vcat(black, white), z=toric3d_slice_checks(lattice, 0))
end

"""
X checks ordered initial white triads, six-edge stars at times `0:T-1`, and
final frontier triads. Every qubit belongs to exactly two checks; the single
global dependency gives rank `(T+2)LyLz-1`.
"""
function star_checks(code::Toric3DEncodedSlice)
    lattice = code.lattice
    m = lattice.Ly * lattice.Lz
    checks = toric3d_encoded_seed_checks(code).x[m+1:end]
    for t in 0:lattice.steps, u in 0:lattice.Ly-1, v in 0:lattice.Lz-1
        support = [data_qubit_index(code, d, t, u, v) for d in 1:3]
        if t < lattice.steps
            append!(support, [data_qubit_index(code, d, t + 1,
                mod(u + shift[1], lattice.Ly), mod(v + shift[2], lattice.Lz))
                for (d, shift) in enumerate(_TORIC3D_SLICE_SHIFTS)])
        end
        push!(checks, support)
    end
    return checks
end

"""Physical plaquettes span the Z checks, including all seed hexagons."""
plaquette_checks(code::Toric3DEncodedSlice) = plaquette_checks(code.lattice)

star_check_matrix(code::Toric3DEncodedSlice; sparse::Bool=false) =
    _check_matrix(star_checks(code), data_qubit_count(code); sparse)
plaquette_check_matrix(code::Toric3DEncodedSlice; sparse::Bool=false) =
    _check_matrix(plaquette_checks(code), data_qubit_count(code); sparse)

"""Transverse X membranes dual to the u and v seed cycles, respectively."""
function toric3d_encoded_logical_x(code::Toric3DEncodedSlice)
    lattice = code.lattice
    return [[data_qubit_index(code, 1, t, 0, v)
             for t in 0:lattice.steps for v in 0:lattice.Lz-1],
            [data_qubit_index(code, 2, t, u, 0)
             for t in 0:lattice.steps for u in 0:lattice.Ly-1]]
end

"""Honeycomb Z cycles winding once around u and v on the initial slice."""
function toric3d_encoded_logical_z(code::Toric3DEncodedSlice)
    lattice = code.lattice
    return [[data_qubit_index(code, d, 0, u, 0)
             for d in (1, 3) for u in 0:lattice.Ly-1],
            [data_qubit_index(code, d, 0, 0, v)
             for d in (2, 3) for v in 0:lattice.Lz-1]]
end

"""The unchanged seven-gate transfer; seed preparation is supplied separately."""
toric3d_encoded_encoder(code::Toric3DEncodedSlice) = toric3d_slice_encoder(code.lattice)

# Row-space membership needs no enumeration of the exponentially large kernel.
function _toric3d_encoded_row_basis(matrix::BitMatrix)
    reduced = copy(matrix)
    pivots = Int[]
    row = 1
    for column in axes(reduced, 2)
        row > size(reduced, 1) && break
        candidate = findfirst(i -> reduced[i, column], row:size(reduced, 1))
        isnothing(candidate) && continue
        source = row + candidate - 1
        reduced[row, :], reduced[source, :] = copy(reduced[source, :]), copy(reduced[row, :])
        for other in row+1:size(reduced, 1)
            reduced[other, column] && (reduced[other, :] .⊻= reduced[row, :])
        end
        push!(pivots, column)
        row += 1
    end
    return (matrix=reduced[1:length(pivots), :], pivots=pivots)
end

function _toric3d_encoded_in_row_space(vector::AbstractVector{Bool}, basis)
    remainder = BitVector(vector)
    for (row, column) in enumerate(basis.pivots)
        remainder[column] && (remainder .⊻= @view basis.matrix[row, :])
    end
    return !any(remainder)
end

function _toric3d_encoded_pulled_to_seed(n, operations, pauli, support,
                                        seed_qubits, seed_z, seed_z_basis)
    x, z = falses(n), falses(n)
    (pauli === :X ? x : z)[support] .= true
    negative = Ref(false)
    for operation in Iterators.reverse(operations)
        operation.gate in (:H, :CNOT) || return false
        _conjugate_pauli!(x, z, negative, operation)
    end
    any(@view x[seed_qubits+1:end]) && return false
    xs, zs = @view(x[1:seed_qubits]), @view(z[1:seed_qubits])
    # In the logical-plus seed, X stabilizers comprise ker(Q), while Z
    # stabilizers comprise row(Q). All fresh qubits remain |0> at the input.
    any(isodd(count(@view(seed_z[row, :]) .& xs)) for row in axes(seed_z, 1)) && return false
    _toric3d_encoded_in_row_space(zs, seed_z_basis) || return false
    overlap = count(xs .& zs)
    iseven(overlap) || return false
    # The Hermitian binary Pauli convention contributes (-1)^(overlap/2)
    # when the commuting pure-X and pure-Z seed stabilizers are multiplied.
    return negative[] == isodd(overlap ÷ 2)
end

"""
Verify every prepared prefix against a honeycomb logical-`|++>` seed and
fresh `|0>` targets using signed binary Clifford conjugation. Checks include
all physical stabilizers, both logical X membranes, and active hexagons.
"""
function verify_toric3d_encoded_encoder(code::Toric3DEncodedSlice,
        encoder::Toric3DSliceEncoder=toric3d_encoded_encoder(code))
    lattice = code.lattice
    actual = encoder.code
    (actual.steps, actual.Ly, actual.Lz) ==
        (lattice.steps, lattice.Ly, lattice.Lz) || return false
    m = lattice.Ly * lattice.Lz
    encoder.checkpoints == [7m*t for t in 1:lattice.steps] || return false
    length(encoder.operations) == encoder.checkpoints[end] || return false
    n = data_qubit_count(code)
    for operation in encoder.operations
        required = operation.gate === :H ? 1 : operation.gate === :CNOT ? 2 : 0
        required > 0 && length(operation.qubits) == required || return false
        all(q -> 1 <= q <= n, operation.qubits) || return false
        length(unique(operation.qubits)) == required || return false
    end
    seed_z = toric3d_slice_check_matrix(lattice)
    basis = _toric3d_encoded_row_basis(seed_z)
    for t in 1:lattice.steps
        prefix = Toric3DEncodedSlice(t, lattice.Ly, lattice.Lz)
        operations = @view encoder.operations[1:encoder.checkpoints[t]]
        all(all(q -> q <= data_qubit_count(prefix), operation.qubits)
            for operation in operations) || return false
        for (pauli, checks) in ((:X, star_checks(prefix)),
                                (:Z, plaquette_checks(prefix)),
                                (:X, toric3d_encoded_logical_x(prefix)),
                                (:Z, toric3d_slice_checks(prefix.lattice, t)))
            all(_toric3d_encoded_pulled_to_seed(n, operations, pauli, support,
                3m, seed_z, basis) for support in checks) || return false
        end
    end
    return true
end
