"""An encoded diagonal prefix and its ordered local-to-full qubit map."""
struct Toric3DEncodedPrefix
    code::Toric3DEncodedSlice
    qubits::Vector{Int}
end

function toric3d_encoded_prefix(code::Toric3DEncodedSlice, t::Integer)
    1 <= t <= code.lattice.steps ||
        throw(ArgumentError("encoded prefix checkpoint is out of range"))
    prefix = Toric3DEncodedSlice(t, code.lattice.Ly, code.lattice.Lz)
    qubits = [data_qubit_index(code,
        data_qubit_coordinate(prefix, q)...) for q in 1:data_qubit_count(prefix)]
    return Toric3DEncodedPrefix(prefix, qubits)
end

"""
The seed or three-sheet star-syndrome window at checkpoint `r`.

`checks` use full-code qubit indices. `qubits` orders the permitted recovery
support, and `matrix` is the corresponding check-by-recovery incidence matrix.
"""
struct Toric3DZFeedbackWindow
    code::Toric3DEncodedSlice
    r::Int
    checks::Vector{Vector{Int}}
    qubits::Vector{Int}
    matrix::BitMatrix
end

function _toric3d_map_prefix_checks(prefix::Toric3DEncodedPrefix,
                                    checks::AbstractVector)
    return [[prefix.qubits[q] for q in check] for check in checks]
end

function _toric3d_restricted_check_matrix(checks::Vector{Vector{Int}},
                                          qubits::Vector{Int})
    columns = Dict(q => column for (column, q) in enumerate(qubits))
    matrix = falses(length(checks), length(qubits))
    for (row, check) in enumerate(checks), q in check
        column = get(columns, q, 0)
        column == 0 || (matrix[row, column] = !matrix[row, column])
    end
    return matrix
end

function toric3d_z_feedback_window(code::Toric3DEncodedSlice, r::Integer)
    0 <= r <= code.lattice.steps ||
        throw(ArgumentError("Z-feedback checkpoint is out of range"))
    checkpoint = Int(r)
    if checkpoint == 0
        checks = deepcopy(toric3d_encoded_seed_checks(code).x)
        qubits = toric3d_slice_qubits(code.lattice, 0)
    else
        prefix = toric3d_encoded_prefix(code, checkpoint)
        prefix_checks = star_checks(prefix.code)
        cells = code.lattice.Ly * code.lattice.Lz
        selected = checkpoint == 1 ? eachindex(prefix_checks) :
                   length(prefix_checks) - 3cells + 1:length(prefix_checks)
        checks = _toric3d_map_prefix_checks(prefix, prefix_checks[selected])
        qubits = vcat(toric3d_slice_qubits(code.lattice, checkpoint - 1),
                      toric3d_slice_qubits(code.lattice, checkpoint))
    end
    matrix = _toric3d_restricted_check_matrix(checks, qubits)
    return Toric3DZFeedbackWindow(code, checkpoint, checks, qubits, matrix)
end

function toric3d_z_window_syndrome(window::Toric3DZFeedbackWindow,
                                    frame::AbstractVector{Bool})
    length(frame) == data_qubit_count(window.code) ||
        throw(DimensionMismatch("frame length does not match the full encoded code"))
    return BitVector(isodd(count(@view frame[check])) for check in window.checks)
end
