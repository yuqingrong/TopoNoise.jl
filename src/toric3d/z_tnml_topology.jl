"""
The p-independent reduced circuit topology used by the Z-only TN decoder.

`raw_detector_matrix` is `D = H*F (mod 2)` and `raw_logical_matrix` is
`B = X*F (mod 2)`, where columns label raw circuit faults. Raw faults with the
same detector endpoints and two-bit logical signature are represented by one
effective edge. `original_columns` preserves the exact raw-fault provenance of
every effective edge.
"""
struct Toric3DZTNTopology
    raw_detector_matrix::SparseArrays.SparseMatrixCSC{UInt8,Int}
    raw_logical_matrix::SparseArrays.SparseMatrixCSC{UInt8,Int}
    detector_endpoints::Vector{NTuple{2,Int}}
    logical_masks::Vector{UInt8}
    multiplicities::Vector{Int}
    original_columns::Vector{Vector{Int}}

    function Toric3DZTNTopology(
            raw_detector_matrix::SparseArrays.SparseMatrixCSC{UInt8,Int},
            raw_logical_matrix::SparseArrays.SparseMatrixCSC{UInt8,Int},
            detector_endpoints::Vector{NTuple{2,Int}},
            logical_masks::Vector{UInt8}, multiplicities::Vector{Int},
            original_columns::Vector{Vector{Int}})
        _validate_toric3d_z_tn_topology(raw_detector_matrix, raw_logical_matrix,
            detector_endpoints, logical_masks, multiplicities, original_columns)
        return new(copy(raw_detector_matrix), copy(raw_logical_matrix),
            copy(detector_endpoints), copy(logical_masks), copy(multiplicities),
            copy.(original_columns))
    end
end

function _validate_toric3d_z_tn_topology(raw_detector_matrix, raw_logical_matrix,
        detector_endpoints, logical_masks, multiplicities, original_columns)
    group_count = length(detector_endpoints)
    all(==(group_count), (length(logical_masks), length(multiplicities),
        length(original_columns))) || throw(DimensionMismatch(
        "effective-edge topology arrays must have equal lengths"))
    length(Set(zip(detector_endpoints, logical_masks))) == group_count ||
        throw(ArgumentError(
            "effective-edge endpoint and logical-mask signatures must be unique"))
    size(raw_logical_matrix, 1) == 2 || throw(DimensionMismatch(
        "the raw logical matrix must have exactly two rows"))
    size(raw_detector_matrix, 2) == size(raw_logical_matrix, 2) ||
        throw(DimensionMismatch(
            "raw detector and logical matrices must have equal column counts"))
    detector_count, raw_fault_count = size(raw_detector_matrix)
    detector_count > 0 && raw_fault_count > 0 ||
        throw(ArgumentError("the raw topology must contain detectors and faults"))
    all(==(0x01), SparseArrays.nonzeros(raw_detector_matrix)) ||
        throw(ArgumentError("the raw detector matrix must be binary"))
    all(==(0x01), SparseArrays.nonzeros(raw_logical_matrix)) ||
        throw(ArgumentError("the raw logical matrix must be binary"))

    seen_columns = falses(raw_fault_count)
    detector_rows = SparseArrays.rowvals(raw_detector_matrix)
    for group in 1:group_count
        endpoint_1, endpoint_2 = detector_endpoints[group]
        1 <= endpoint_1 < endpoint_2 <= detector_count || throw(ArgumentError(
            "effective-edge endpoints must be ordered, distinct, and in bounds"))
        logical_masks[group] <= 0x03 ||
            throw(ArgumentError("logical masks must lie in 0:3"))
        multiplicities[group] > 0 ||
            throw(ArgumentError("effective-edge multiplicities must be positive"))
        multiplicities[group] == length(original_columns[group]) ||
            throw(ArgumentError(
                "each multiplicity must equal its raw-column provenance count"))
        for column in original_columns[group]
            1 <= column <= raw_fault_count ||
                throw(ArgumentError("raw-column provenance is out of bounds"))
            seen_columns[column] &&
                throw(ArgumentError("raw-column provenance must not overlap"))
            seen_columns[column] = true
            rows = detector_rows[SparseArrays.nzrange(raw_detector_matrix, column)]
            length(rows) == 2 && (rows[1], rows[2]) == detector_endpoints[group] ||
                throw(ArgumentError(
                    "raw detector column does not match its effective edge"))
            raw_mask = raw_logical_matrix[1, column] |
                       (raw_logical_matrix[2, column] << 1)
            raw_mask == logical_masks[group] || throw(ArgumentError(
                "raw logical column does not match its effective edge"))
        end
    end
    all(seen_columns) ||
        throw(ArgumentError("raw-column provenance must be a complete partition"))
    return nothing
end

"""
    build_toric3d_z_tn_topology(model)

Compile the exact propagated circuit-fault matrix into a reduced detector
graph. Each raw fault must have exactly two detector endpoints. Parallel raw
faults are combined only when both endpoints and both logical bits agree.
"""
function build_toric3d_z_tn_topology(model::Toric3DZFaultModel)
    code = model.code
    detector_checks = UInt8.(star_check_matrix(code; sparse=true))
    logical_checks = UInt8.(_check_matrix(toric3d_encoded_logical_x(code),
        data_qubit_count(code); sparse=true))
    raw_detector_matrix = SparseArrays.sparse(
        UInt8.(isodd.(detector_checks * model.matrix)))
    raw_logical_matrix = SparseArrays.sparse(
        UInt8.(isodd.(logical_checks * model.matrix)))

    detector_endpoints = NTuple{2,Int}[]
    logical_masks = UInt8[]
    original_columns = Vector{Int}[]
    group_indices = Dict{Tuple{NTuple{2,Int},UInt8},Int}()
    row_indices = SparseArrays.rowvals(raw_detector_matrix)
    for column in axes(raw_detector_matrix, 2)
        rows = row_indices[SparseArrays.nzrange(raw_detector_matrix, column)]
        length(rows) == 2 || error(
            "raw Z fault column $column must have exactly two detector endpoints")
        endpoints = (rows[1], rows[2])
        logical_mask = raw_logical_matrix[1, column] |
                       (raw_logical_matrix[2, column] << 1)
        signature = (endpoints, logical_mask)
        group = get(group_indices, signature, 0)
        if iszero(group)
            push!(detector_endpoints, endpoints)
            push!(logical_masks, logical_mask)
            push!(original_columns, [column])
            group_indices[signature] = length(original_columns)
        else
            push!(original_columns[group], column)
        end
    end
    multiplicities = length.(original_columns)
    return Toric3DZTNTopology(raw_detector_matrix, raw_logical_matrix,
        detector_endpoints, logical_masks, multiplicities, original_columns)
end

"""
    toric3d_z_tn_edge_probabilities(topology, p)

Return the probability that each effective edge has odd parity. An edge with
`k` identical raw faults has `q_k = (1 - (1 - 2p)^k)/2`.
"""
function toric3d_z_tn_edge_probabilities(topology::Toric3DZTNTopology, p::Real)
    _validate_fault_rate(p, "p")
    rate = Float64(p)
    return [-expm1(multiplicity * log1p(-2rate)) / 2
            for multiplicity in topology.multiplicities]
end

