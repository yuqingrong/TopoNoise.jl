"""Correction bonds output by the MWPM decoder.

`horizontal` and `vertical` share the shapes of `VirtualBondErrors.horizontal_internal`
`(rows, cols - 1)` and `.vertical_internal` `(rows - 1, cols)`. Each entry is `true`
if the decoder chose to flip the corresponding internal virtual bond as part of the
correction chain.
"""
struct Correction
    horizontal::BitMatrix
    vertical::BitMatrix
end

Correction(rows::Integer, cols::Integer) =
    Correction(falses(Int(rows), Int(cols) - 1), falses(Int(rows) - 1, Int(cols)))

residual_errors(errors::DataEdgeErrors, correction::Correction) =
    _residual_errors(errors, correction)

function logical_failure(
        errors::DataEdgeErrors, correction::Correction;
        sector::Symbol=:north_south)
    return _data_edge_logical_failure(errors, correction; sector=sector)
end

"""
    plaquette_syndrome(mismatches) -> BitMatrix

Return the `(R-1) x (C-1)` interior plaquette-parity syndrome of a mismatch
record. Entry `[r, c]` is `true` iff an odd number of the four bonds
surrounding plaquette `(r, c)` are mismatched.
"""
function plaquette_syndrome(mismatches)
    horizontal = mismatches.horizontal
    vertical = mismatches.vertical
    rows = size(horizontal, 1)
    cols = size(vertical, 2)
    size(horizontal) == (rows, cols - 1) || throw(DimensionMismatch(
        "horizontal mismatches must have size ($rows, $(cols - 1))"))
    size(vertical) == (rows - 1, cols) || throw(DimensionMismatch(
        "vertical mismatches must have size ($(rows - 1), $cols)"))
    syndrome = BitMatrix(undef, rows - 1, cols - 1)
    for r in 1:(rows - 1), c in 1:(cols - 1)
        syndrome[r, c] = xor(
            horizontal[r, c], horizontal[r + 1, c],
            vertical[r, c], vertical[r, c + 1])
    end
    return syndrome
end

const _BOUNDARY_MODES = (:north_south, :east_west)

_validate_boundary(boundary::Symbol) = (
    boundary in _BOUNDARY_MODES || throw(ArgumentError(
        "boundary must be :north_south or :east_west, got $boundary"));
    boundary)

const _pymatching_ref = Ref{Py}()
const _matching_cache = Dict{Tuple{Int,Int,Symbol},Py}()

_pymatching() = (isassigned(_pymatching_ref) ||
    (_pymatching_ref[] = pyimport("pymatching")); _pymatching_ref[])

_h_qubit(r::Int, c::Int, cols::Int) = (r - 1) * (cols - 1) + (c - 1)
_v_qubit(r::Int, c::Int, rows::Int, cols::Int) =
    rows * (cols - 1) + (r - 1) * cols + (c - 1)

_num_qubits(rows::Int, cols::Int) = rows * (cols - 1) + (rows - 1) * cols

pyset(v::AbstractVector) = pybuiltins.set(pylist(collect(v)))
pyset(s::AbstractSet) = pyset(collect(s))
pyint(n::Integer) = pybuiltins.int(n)

# Build the PyMatching graph for the surface code:
#   :north_south — N and S rough (bonds H[1,c] and H[R,c] escape to boundary);
#                  E and W smooth (bonds V[r,1] and V[r,C] escape to boundary).
#                  Logical X = column of V bonds at c=c_ref → detector = a row
#                  cut so fault_ids of V[r_ref, c] carry logical.
#   :east_west   — transposed.
# Each internal bond also gets a unique fault_id equal to its qubit index so
# the decoder's correction can be recovered per-bond.
function _build_matching(rows::Int, cols::Int, boundary::Symbol)
    pm = _pymatching()
    R, C = rows, cols
    matching = pm.Matching()
    n_qubits = _num_qubits(R, C)
    logical_id = n_qubits  # fault id for the logical observable
    interior_check(r, c) = (r - 1) * (C - 1) + (c - 1)  # 0-based

    r_cut_ns = boundary === :north_south ? max(1, cld(R - 1, 2)) : 0
    c_cut_ew = boundary === :east_west   ? max(1, cld(C - 1, 2)) : 0

    # Horizontal bonds H[r, c]: r in 1..R, c in 1..C-1.
    # H[r,c] connects plaquettes (r-1, c) and (r, c). Missing when r==1 (north)
    # or r==R (south).
    for r in 1:R, c in 1:(C - 1)
        qubit = _h_qubit(r, c, C)
        checks = Int[]
        r > 1 && push!(checks, interior_check(r - 1, c))
        r < R && push!(checks, interior_check(r, c))
        fault_ids = [qubit]
        boundary === :east_west && c == c_cut_ew && push!(fault_ids, logical_id)
        if length(checks) == 2
            matching.add_edge(pyint(checks[1]), pyint(checks[2]),
                              fault_ids=pyset(fault_ids), weight=1.0,
                              merge_strategy="independent")
        else
            matching.add_boundary_edge(pyint(checks[1]),
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="independent")
        end
    end

    # Vertical bonds V[r, c]: r in 1..R-1, c in 1..C.
    # V[r,c] connects plaquettes (r, c-1) and (r, c). Missing when c==1 (west)
    # or c==C (east).
    for r in 1:(R - 1), c in 1:C
        qubit = _v_qubit(r, c, R, C)
        checks = Int[]
        c > 1 && push!(checks, interior_check(r, c - 1))
        c < C && push!(checks, interior_check(r, c))
        fault_ids = [qubit]
        boundary === :north_south && r == r_cut_ns &&
            push!(fault_ids, logical_id)
        if length(checks) == 2
            matching.add_edge(pyint(checks[1]), pyint(checks[2]),
                              fault_ids=pyset(fault_ids), weight=1.0,
                              merge_strategy="independent")
        else
            matching.add_boundary_edge(pyint(checks[1]),
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="independent")
        end
    end
    return matching
end

function _matching(rows::Int, cols::Int, boundary::Symbol)
    _validate_boundary(boundary)
    key = (rows, cols, boundary)
    haskey(_matching_cache, key) && return _matching_cache[key]
    _matching_cache[key] = _build_matching(rows, cols, boundary)
    return _matching_cache[key]
end

const _open_matching_cache = Dict{Tuple{Int,Int,Symbol},Py}()

function _validate_open_syndrome(
        model::OpenCodeCapacityModel, syndrome::BitMatrix)
    expected = (model.rows - 1, model.cols - 1)
    size(syndrome) == expected || throw(DimensionMismatch(
        "syndrome must have size $expected"))
    return syndrome
end

function _syndrome_vector(syndrome::BitMatrix)
    rows, cols = size(syndrome)
    flattened = Vector{Int8}(undef, rows * cols)
    for row in 1:rows, col in 1:cols
        flattened[(row - 1) * cols + col] = syndrome[row, col] ? 1 : 0
    end
    return pyimport("numpy").asarray(flattened)
end

function _build_open_code_matching(
        model::OpenCodeCapacityModel, sector::Symbol)
    _validate_boundary(sector)
    rows, cols = model.rows, model.cols
    matching = _pymatching().Matching()
    logical_id = _num_qubits(rows, cols)
    cut = logical_cut(model; sector=sector)
    boundary_nodes = Int[]
    next_boundary_node = (rows - 1) * (cols - 1)
    interior_check(row, col) = (row - 1) * (cols - 1) + (col - 1)

    for row in 1:rows, col in 1:(cols - 1)
        qubit = _h_qubit(row, col, cols)
        checks = Int[]
        row > 1 && push!(checks, interior_check(row - 1, col))
        row < rows && push!(checks, interior_check(row, col))
        fault_ids = cut.horizontal[row, col] ? [qubit, logical_id] : [qubit]
        if length(checks) == 2
            matching.add_edge(pyint(checks[1]), pyint(checks[2]);
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="disallow")
        else
            push!(boundary_nodes, next_boundary_node)
            matching.add_edge(pyint(only(checks)), pyint(next_boundary_node);
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="disallow")
            next_boundary_node += 1
        end
    end

    for row in 1:(rows - 1), col in 1:cols
        qubit = _v_qubit(row, col, rows, cols)
        checks = Int[]
        col > 1 && push!(checks, interior_check(row, col - 1))
        col < cols && push!(checks, interior_check(row, col))
        fault_ids = cut.vertical[row, col] ? [qubit, logical_id] : [qubit]
        if length(checks) == 2
            matching.add_edge(pyint(checks[1]), pyint(checks[2]);
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="disallow")
        else
            push!(boundary_nodes, next_boundary_node)
            matching.add_edge(pyint(only(checks)), pyint(next_boundary_node);
                fault_ids=pyset(fault_ids), weight=1.0,
                merge_strategy="disallow")
            next_boundary_node += 1
        end
    end
    matching.set_boundary_nodes(pyset(boundary_nodes))
    return matching
end

function _open_code_matching(
        model::OpenCodeCapacityModel, sector::Symbol;
        error_rate::Float64)
    _validate_boundary(sector)
    key = (model.rows, model.cols, sector)
    return get!(_open_matching_cache, key) do
        _build_open_code_matching(model, sector)
    end
end

function _correction_from_faults(
        model::OpenCodeCapacityModel, predicted)::Correction
    prediction = pyconvert(Vector{Int}, predicted)
    required = _num_qubits(model.rows, model.cols)
    length(prediction) >= required || throw(ErrorException(
        "pymatching returned $(length(prediction)) fault predictions; expected " *
        "at least $required"))
    correction = Correction(model.rows, model.cols)
    for row in 1:model.rows, col in 1:(model.cols - 1)
        correction.horizontal[row, col] =
            prediction[_h_qubit(row, col, model.cols) + 1] == 1
    end
    for row in 1:(model.rows - 1), col in 1:model.cols
        correction.vertical[row, col] =
            prediction[_v_qubit(row, col, model.rows, model.cols) + 1] == 1
    end
    return correction
end

"""Decode an open-patch plaquette syndrome without access to true errors."""
function decode_syndrome(
        model::OpenCodeCapacityModel, syndrome::BitMatrix;
        sector::Symbol=:north_south, error_rate::Real=0.1)::Correction
    _validate_open_syndrome(model, syndrome)
    _validate_boundary(sector)
    0 <= error_rate < 0.5 ||
        throw(ArgumentError("error_rate must be in [0, 0.5)"))
    any(syndrome) || return Correction(model.rows, model.cols)
    error_rate > 0 ||
        throw(ArgumentError("nonzero syndrome is impossible at p=0"))
    matching = _open_code_matching(
        model, sector; error_rate=Float64(error_rate))
    predicted = matching.decode(_syndrome_vector(syndrome))
    return _correction_from_faults(model, predicted)
end

"""
    decode_uf(model, mismatches; boundary=:north_south) -> Correction

Decode a mismatch record with PyMatching's blossom-based MWPM. `:north_south`
treats N/S as rough (bond-column endpoints on H[1,:] and H[R,:] escape to
the boundary) and E/W as smooth (V[:,1] and V[:,C] also escape to the
boundary — smooth-boundary chains are freely satisfiable and produce no
detectable defect). Logical X operator = column of V bonds; logical Z
detector = parity of V[r_ref, :] with r_ref at the mid-row.
`:east_west` is the transposed geometry.
"""
function decode_uf(
        model::ToricCodeTrajectoryModel, mismatches;
        boundary::Symbol=:north_south)
    _validate_boundary(boundary)
    rows, cols = model.rows, model.cols
    correction = Correction(rows, cols)
    (rows < 2 || cols < 2) && return correction

    R, C = rows, cols
    interior = plaquette_syndrome(mismatches)
    syndrome_vec = zeros(Int8, (R - 1) * (C - 1))
    for r in 1:(R - 1), c in 1:(C - 1)
        syndrome_vec[(r - 1) * (C - 1) + c] = interior[r, c] ? 1 : 0
    end
    any(syndrome_vec .!= 0) || return correction

    matching = _matching(R, C, boundary)
    np = pyimport("numpy")
    syndrome_py = np.asarray(syndrome_vec)
    corr_py = matching.decode(syndrome_py)
    corr_vec = pyconvert(Vector{Int}, corr_py)
    n_qubits = _num_qubits(R, C)
    length(corr_vec) >= n_qubits || throw(ErrorException(
        "pymatching returned $(length(corr_vec)) fault predictions; expected " *
        "at least $n_qubits"))

    for r in 1:R, c in 1:(C - 1)
        correction.horizontal[r, c] = corr_vec[_h_qubit(r, c, C) + 1] == 1
    end
    for r in 1:(R - 1), c in 1:C
        correction.vertical[r, c] = corr_vec[_v_qubit(r, c, R, C) + 1] == 1
    end
    return correction
end

"""
    logical_failure(errors_true, correction; boundary=:north_south) -> Bool

Parity of internal V-bond flips (error XOR correction) crossing a fixed
reference row for `:north_south`, or of H-bond flips crossing a fixed
reference column for `:east_west`.
"""
function logical_failure(
        errors_true::VirtualBondErrors, correction::Correction;
        boundary::Symbol=:north_south)
    _validate_boundary(boundary)
    if boundary === :north_south
        rows_int, cols_int = size(errors_true.vertical_internal)
        size(correction.vertical) == (rows_int, cols_int) ||
            throw(DimensionMismatch(
                "correction.vertical must match vertical_internal shape"))
        r_ref = max(1, cld(rows_int, 2))
        parity = false
        for c in 1:cols_int
            parity = xor(parity,
                errors_true.vertical_internal[r_ref, c],
                correction.vertical[r_ref, c])
        end
        return parity
    else
        rows_int, cols_int = size(errors_true.horizontal_internal)
        size(correction.horizontal) == (rows_int, cols_int) ||
            throw(DimensionMismatch(
                "correction.horizontal must match horizontal_internal shape"))
        c_ref = max(1, cld(cols_int, 2))
        parity = false
        for r in 1:rows_int
            parity = xor(parity,
                errors_true.horizontal_internal[r, c_ref],
                correction.horizontal[r, c_ref])
        end
        return parity
    end
end

"""One decoded trajectory: the decoder correction plus its logical verdict."""
struct DecodedTrajectory
    correction::Correction
    boundary::Symbol
    logical_failure::Bool
end

"""
    decode_trajectory(trajectory; boundary=:north_south) -> DecodedTrajectory
"""
function decode_trajectory(
        trajectory::ToricCodeTrajectory; boundary::Symbol=:north_south)
    rows, cols, _ = size(trajectory.measurements)
    model = ToricCodeTrajectoryModel(rows, cols)
    mismatches = bond_mismatches(trajectory)
    correction = decode_uf(model, mismatches; boundary=boundary)
    failed = logical_failure(
        trajectory.errors, correction; boundary=boundary)
    return DecodedTrajectory(correction, boundary, failed)
end
