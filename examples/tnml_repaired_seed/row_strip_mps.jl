if !isdefined(@__MODULE__, :TNMLRowStripPlan)
    include(joinpath(@__DIR__, "row_strip_plan.jl"))
end

"""
Experimental row-strip boundary-MPS decoder.

The global state after each strip is an MPS over the active detector parities,
in the physical column order compiled by `TNMLRowStripPlan`, followed by one
four-state logical-class site.  A row map introduces all independent edges
owned by that row, fixes the observed syndrome on the retired detector row,
and emits only the detector parities that remain active.

The dense row map remains as an exact-small reference.  One open-chain
finite-state row MPO carries incoming cross-column edge bits, two logical
bits, the explicit periodic-seam assignment, and one accumulated
source-parity bit.

The lower-bond reference uses four logical Walsh/Fourier sectors.  A sector
MPS is a signed partial partition function over active detector assignments,
not a probability distribution.  The logical carrier disappears; after
fixing the periodic seam, the row automaton carries only the incoming cross
bits and source parity, with maximum dimension eight for this stencil.  Each
signed sector is therefore normalized by a positive Euclidean MPS norm rather
than by its sum of amplitudes.
"""
module TNMLRowStripMPS

using ITensors, ITensorMPS, TopoNoise
using ..TNMLRowStripPlan

export DenseRowMap, RowMPO, FourierRowMPO, RowStripMPSDecoder,
       compile_dense_row_map,
       compile_row_mpo, compile_fixed_seam_row_mpos, row_mpo_action,
       compile_fourier_fixed_seam_row_mpos,
       build_row_strip_mps_decoder, decode_row_strip_mps,
       decode_row_strip_mps!

const MAX_DENSE_MAP_ENTRIES = 16_000_000
const MAX_DENSE_DETECTORS = Sys.WORD_SIZE - 3

struct DenseRowMap
    strip_id::Int
    input_detectors::Vector{Int}
    output_detectors::Vector{Int}
    matrix::Matrix{Float64}
    transition_count::Int
end

struct RowMPO
    strip_id::Int
    union_detectors::Vector{Int}
    operator::ITensorMPS.MPO
    cross_offsets::Vector{NTuple{3,Int}}
    maximum_automaton_dimension::Int
    seam::Union{Nothing,Int}
end

struct FourierRowMPO
    strip_id::Int
    union_detectors::Vector{Int}
    operator::ITensorMPS.MPO
    cross_offsets::Vector{NTuple{3,Int}}
    maximum_automaton_dimension::Int
    seam::Int
    sector::Int
end

struct RowStripMPSDecoder
    plan::RowStripPlan
    probabilities::Vector{Float64}
    detector_sites::Dict{Int,Index}
    logical_site::Index
    scalar_site::Index
    maxdim::Int
    cutoff::Float64
    exact_bond_guard::Int
end

struct _BoundaryMPS
    state::ITensorMPS.MPS
    detectors::Vector{Int}
end

struct _SectorBoundaryMPS
    state::ITensorMPS.MPS
    detectors::Vector{Int}
end

function _validate_probabilities(plan, probabilities)
    length(probabilities) == length(plan.edges) || throw(DimensionMismatch(
        "incorrect effective-edge probability count"))
    values = Float64.(probabilities)
    all(value -> isfinite(value) && 0 <= value <= 1, values) ||
        throw(ArgumentError("edge probabilities must be finite and in [0,1]"))
    return values
end

function _validate_syndrome(plan, syndrome; require_even::Bool=true)
    length(syndrome) == length(plan.detector_coordinates) ||
        throw(DimensionMismatch("incorrect syndrome length"))
    all(bit -> bit == 0 || bit == 1, syndrome) ||
        throw(ArgumentError("syndrome must be binary"))
    !require_even || iseven(count(!iszero, syndrome)) ||
        throw(ArgumentError("odd parity syndrome is impossible"))
    return BitVector(syndrome)
end

function _physical_detector_order(plan, detectors)
    sort!(collect(detectors); by=detector -> begin
        coordinate = plan.detector_coordinates[detector]
        (coordinate.column, coordinate.id)
    end)
end

function _workspace(plan, strip)
    detectors = Int[]
    append!(detectors, (slot.detector for slot in strip.physical_before))
    append!(detectors, strip.detectors)
    for edge_id in strip.owned_edges
        append!(detectors, plan.edges[edge_id].endpoints)
    end
    return _physical_detector_order(plan, unique(detectors))
end

function _add_weight!(distribution, key, weight)
    iszero(weight) && return distribution
    distribution[key] = get(distribution, key, 0.0) + weight
    return distribution
end

function _edge_delta_distribution(plan, strip, probabilities, positions)
    distribution = Dict{Tuple{UInt64,UInt8},Float64}((0x00, 0x00) => 1.0)
    for edge_id in strip.owned_edges
        edge = plan.edges[edge_id]
        probability = probabilities[edge_id]
        delta = (UInt64(1) << (positions[edge.endpoints[1]] - 1)) |
                (UInt64(1) << (positions[edge.endpoints[2]] - 1))
        updated = Dict{Tuple{UInt64,UInt8},Float64}()
        for ((mask, logical), weight) in distribution
            _add_weight!(updated, (mask, logical), (1 - probability) * weight)
            _add_weight!(updated, (xor(mask, delta),
                xor(logical, edge.logical_mask)), probability * weight)
        end
        distribution = updated
    end
    return distribution
end

function _scatter_masks(detectors, positions)
    count = length(detectors)
    mapped = zeros(UInt64, 1 << count)
    deltas = [UInt64(1) << (positions[detector] - 1)
              for detector in detectors]
    for mask in 1:length(mapped)-1
        bit = trailing_zeros(mask) + 1
        mapped[mask + 1] = xor(mapped[xor(mask, 1 << (bit - 1)) + 1],
            deltas[bit])
    end
    return mapped
end

function _gather_mask(mask::UInt64, detectors, positions)
    output = 0
    for (bit, detector) in enumerate(detectors)
        iszero(mask & (UInt64(1) << (positions[detector] - 1))) ||
            (output |= 1 << (bit - 1))
    end
    return output
end

function _dense_dimension(detector_count::Integer)
    0 <= detector_count <= MAX_DENSE_DETECTORS || throw(ArgumentError(
        "dense semantic prototype cannot index $detector_count detector bits"))
    return 4 * (1 << detector_count)
end

function _check_dense_map_size(input_count, output_count)
    input_dimension = _dense_dimension(input_count)
    output_dimension = _dense_dimension(output_count)
    big(input_dimension) * output_dimension <= MAX_DENSE_MAP_ENTRIES ||
        throw(ArgumentError(
            "dense semantic prototype exceeds $(MAX_DENSE_MAP_ENTRIES) map entries"))
    return input_dimension, output_dimension
end

"""
    compile_dense_row_map(plan, strip_id, probabilities, syndrome)

Materialize the exact conditioned transfer for one row.  Input and output
bases use detector bits as low-order bits in `strip.physical_*` order and the
two-bit logical class above them.  The map includes periodic seam edges by
their compiled detector IDs; there is no separate boundary convention to
guess in the numerical code.
"""
function compile_dense_row_map(plan::RowStripPlan, strip_id::Integer,
        probabilities, syndrome)
    1 <= strip_id <= length(plan.strips) ||
        throw(BoundsError(plan.strips, strip_id))
    values = _validate_probabilities(plan, probabilities)
    observed = _validate_syndrome(plan, syndrome; require_even=false)
    strip = plan.strips[Int(strip_id)]
    input_detectors = [slot.detector for slot in strip.physical_before]
    output_detectors = [slot.detector for slot in strip.physical_after]
    input_dimension, output_dimension = _check_dense_map_size(
        length(input_detectors), length(output_detectors))

    workspace = _workspace(plan, strip)
    length(workspace) <= 63 || throw(ArgumentError(
        "dense semantic prototype supports at most 63 workspace detectors"))
    positions = Dict(detector => position
        for (position, detector) in enumerate(workspace))
    input_masks = _scatter_masks(input_detectors, positions)
    transitions = _edge_delta_distribution(plan, strip, values, positions)
    matrix = zeros(Float64, output_dimension, input_dimension)

    for input_logical in 0:3, input_mask in 0:(1 << length(input_detectors)) - 1
        initial = input_masks[input_mask + 1]
        column = input_mask + (input_logical << length(input_detectors)) + 1
        for ((delta, logical_delta), weight) in transitions
            combined = xor(initial, delta)
            all(detector ->
                !iszero(combined &
                    (UInt64(1) << (positions[detector] - 1))) ==
                    observed[detector], strip.retired_detectors) || continue
            output_mask = _gather_mask(
                combined, output_detectors, positions)
            output_logical = xor(input_logical, Int(logical_delta))
            row = output_mask +
                  (output_logical << length(output_detectors)) + 1
            matrix[row, column] += weight
        end
    end
    return DenseRowMap(Int(strip_id), input_detectors, output_detectors,
        matrix, length(transitions))
end

function _maximum_exact_bond(detector_count)
    dimensions = vcat(fill(2, detector_count), 4)
    total = prod(big, dimensions)
    left = big(1)
    maximum_bond = big(1)
    for position in 1:length(dimensions)-1
        left *= dimensions[position]
        maximum_bond = max(maximum_bond, min(left, total ÷ left))
    end
    maximum_bond <= typemax(Int) || return typemax(Int)
    return Int(maximum_bond)
end

"""
    build_row_strip_mps_decoder(topology, L, p; maxdim=4096, cutoff=0)

Build a row-strip decoder.  The dense reference remains guarded inside
`compile_dense_row_map`; decoder construction itself supports the open-chain
row-MPO path at larger L.
"""
function build_row_strip_mps_decoder(topology, L::Integer, p::Real;
        maxdim::Integer=4096, cutoff::Real=0.0)
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    isfinite(cutoff) && 0 <= cutoff < 1 ||
        throw(ArgumentError("cutoff must be finite and in [0,1)"))
    plan = compile_repaired_seed_transfer_plan(topology, L)
    probabilities = toric3d_z_tn_edge_probabilities(topology, p)
    detector_sites = Dict(coordinate.id => Index(2,
        "Detector,id=$(coordinate.id),plane=$(coordinate.plane)," *
        "row=$(coordinate.row),column=$(coordinate.column)")
        for coordinate in plan.detector_coordinates)
    logical_site = Index(4, "LogicalClass")
    scalar_site = Index(1, "FourierSectorScalar")
    exact_bond_guard = maximum(
        (_maximum_exact_bond(length(strip.physical_after))
         for strip in plan.strips); init=1)
    return RowStripMPSDecoder(plan, probabilities, detector_sites,
        logical_site, scalar_site, Int(maxdim), Float64(cutoff),
        exact_bond_guard)
end

_boundary_state_count(cross_count) = begin
    assignments = 1 << cross_count
    4 * assignments * assignments
end

_work_state_count(cross_count) = 2 * _boundary_state_count(cross_count)

_fixed_boundary_state_count(cross_count) = 4 * (1 << cross_count)
_fixed_work_state_count(cross_count) = 2 * _fixed_boundary_state_count(cross_count)

_row_boundary_state_count(cross_count, seam) = isnothing(seam) ?
    _boundary_state_count(cross_count) : _fixed_boundary_state_count(cross_count)
_row_work_state_count(cross_count, seam) = isnothing(seam) ?
    _work_state_count(cross_count) : _fixed_work_state_count(cross_count)

function _boundary_state_index(incoming, logical, seam, cross_count)
    assignments = 1 << cross_count
    return 1 + incoming + assignments * (logical + 4 * seam)
end

function _decode_boundary_state(index, cross_count)
    assignments = 1 << cross_count
    value = index - 1
    incoming = value % assignments
    value ÷= assignments
    logical = value % 4
    seam = value ÷ 4
    return incoming, logical, seam
end

function _work_state_index(incoming, logical, seam, source, cross_count)
    return 1 + source + 2 *
        (_boundary_state_index(incoming, logical, seam, cross_count) - 1)
end

function _decode_work_state(index, cross_count)
    value = index - 1
    source = value % 2
    incoming, logical, seam = _decode_boundary_state(value ÷ 2 + 1,
        cross_count)
    return incoming, logical, seam, source
end

function _row_boundary_state_index(incoming, logical, carried_seam,
        cross_count, fixed_seam)
    if isnothing(fixed_seam)
        return _boundary_state_index(
            incoming, logical, carried_seam, cross_count)
    end
    carried_seam == fixed_seam || error("fixed seam changed inside row MPO")
    assignments = 1 << cross_count
    return 1 + incoming + assignments * logical
end

function _decode_row_boundary_state(index, cross_count, fixed_seam)
    isnothing(fixed_seam) && return _decode_boundary_state(index, cross_count)
    assignments = 1 << cross_count
    value = index - 1
    incoming = value % assignments
    logical = value ÷ assignments
    return incoming, logical, fixed_seam
end

function _row_work_state_index(incoming, logical, carried_seam, source,
        cross_count, fixed_seam)
    return 1 + source + 2 *
        (_row_boundary_state_index(incoming, logical, carried_seam,
            cross_count, fixed_seam) - 1)
end

function _decode_row_work_state(index, cross_count, fixed_seam)
    value = index - 1
    source = value % 2
    incoming, logical, carried_seam = _decode_row_boundary_state(
        value ÷ 2 + 1, cross_count, fixed_seam)
    return incoming, logical, carried_seam, source
end


function _edge_tables(plan, strip)
    size = plan.L
    sources = Dict(coordinate.column => coordinate.id
        for coordinate in plan.detector_coordinates[strip.detectors])
    length(sources) == size || error("strip source row is incomplete")
    by_column = [Dict{NTuple{3,Int},Int}() for _ in 1:size]
    for edge_id in strip.owned_edges
        edge = plan.edges[edge_id]
        source = plan.detector_coordinates[edge.endpoints[1]]
        haskey(by_column[source.column + 1], edge.offset) && error(
            "row MPO requires one effective edge per source and offset")
        by_column[source.column + 1][edge.offset] = edge_id
    end
    offset_sets = [Set(keys(table)) for table in by_column]
    all(==(offset_sets[1]), offset_sets) || error(
        "row MPO requires the same geometric edge stencil in every column")
    cross_offsets = sort!(collect(filter(offset -> offset[3] == 1,
        offset_sets[1])))
    same_offsets = sort!(collect(filter(offset -> offset[3] == 0,
        offset_sets[1])))
    return sources, by_column, same_offsets, cross_offsets
end

function _edge_choices(edge_ids, plan, probabilities)
    choices = Tuple{Int,Int,Float64}[(0, 0, 1.0)]
    for edge_id in edge_ids
        edge = plan.edges[edge_id]
        probability = probabilities[edge_id]
        updated = Tuple{Int,Int,Float64}[]
        for (parity, logical, weight) in choices
            iszero((1 - probability) * weight) ||
                push!(updated, (parity, logical,
                    (1 - probability) * weight))
            iszero(probability * weight) ||
                push!(updated, (xor(parity, 1),
                    xor(logical, Int(edge.logical_mask)),
                    probability * weight))
        end
        choices = updated
    end
    return choices
end

function _cross_choices(edge_ids, plan, probabilities)
    choices = Tuple{Int,Int,Int,Float64}[(0, 0, 0, 1.0)]
    for (bit, edge_id) in enumerate(edge_ids)
        edge = plan.edges[edge_id]
        probability = probabilities[edge_id]
        updated = Tuple{Int,Int,Int,Float64}[]
        for (mask, parity, logical, weight) in choices
            iszero((1 - probability) * weight) ||
                push!(updated, (mask, parity, logical,
                    (1 - probability) * weight))
            iszero(probability * weight) ||
                push!(updated, (mask | (1 << (bit - 1)), xor(parity, 1),
                    xor(logical, Int(edge.logical_mask)),
                    probability * weight))
        end
        choices = updated
    end
    return choices
end

function _same_edges_by_target(plan, table, same_offsets)
    result = Dict{Int,Vector{Int}}()
    for offset in same_offsets
        edge_id = table[offset]
        destination = plan.edges[edge_id].endpoints[2]
        push!(get!(result, destination, Int[]), edge_id)
    end
    return result
end

function _incoming_types_by_target(plan, previous_table, cross_offsets)
    result = Dict{Int,Vector{Int}}()
    for (bit, offset) in enumerate(cross_offsets)
        edge_id = previous_table[offset]
        destination = plan.edges[edge_id].endpoints[2]
        push!(get!(result, destination, Int[]), bit)
    end
    return result
end

function _incoming_toggle(incoming, bit_positions)
    parity = 0
    for bit in bit_positions
        parity = xor(parity, (incoming >> (bit - 1)) & 1)
    end
    return parity
end

function _add_tensor_entry!(array, indices::Tuple, value)
    iszero(value) || (array[indices...] += value)
    return array
end

function _source_tensor(site, left_link, right_link, cross_count;
        initial::Bool=false, fixed_seam::Union{Nothing,Int}=nothing)
    output = prime(site)
    work_dimension = _row_work_state_count(cross_count, fixed_seam)
    if initial
        array = zeros(Float64, work_dimension, 2, 2)
        seams = isnothing(fixed_seam) ?
            (0:(1 << cross_count)-1) : (fixed_seam:fixed_seam)
        for seam in seams, input in 0:1
            right = _row_work_state_index(
                seam, 0, seam, input, cross_count, fixed_seam)
            array[right, 1, input + 1] += 1.0
        end
        return ITensor(array, right_link, output, site)
    end
    boundary_dimension = _row_boundary_state_count(cross_count, fixed_seam)
    array = zeros(Float64, boundary_dimension, work_dimension, 2, 2)
    for left in 1:boundary_dimension, input in 0:1
        incoming, logical, seam = _decode_row_boundary_state(
            left, cross_count, fixed_seam)
        right = _row_work_state_index(
            incoming, logical, seam, input, cross_count, fixed_seam)
        array[left, right, 1, input + 1] += 1.0
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _working_tensor(site, left_link, right_link, detector,
        same_edges, incoming_types, plan, probabilities, cross_count;
        fixed_seam::Union{Nothing,Int}=nothing)
    output = prime(site)
    work_dimension = _row_work_state_count(cross_count, fixed_seam)
    array = zeros(Float64, work_dimension, work_dimension, 2, 2)
    choices = _edge_choices(same_edges, plan, probabilities)
    for left in 1:work_dimension, input in 0:1
        incoming, logical, seam, source =
            _decode_row_work_state(left, cross_count, fixed_seam)
        incoming_parity = _incoming_toggle(incoming, incoming_types)
        for (same_parity, logical_delta, weight) in choices
            output_bit = xor(input, xor(incoming_parity, same_parity))
            right = _row_work_state_index(incoming,
                xor(logical, logical_delta), seam,
                xor(source, same_parity), cross_count, fixed_seam)
            _add_tensor_entry!(array,
                (left, right, output_bit + 1, input + 1), weight)
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _final_working_tensor(site, left_link, right_link, detector,
        same_edges, incoming_types, cross_edges, source_syndrome,
        plan, probabilities, cross_count;
        fixed_seam::Union{Nothing,Int}=nothing)
    output = prime(site)
    work_dimension = _row_work_state_count(cross_count, fixed_seam)
    boundary_dimension = _row_boundary_state_count(cross_count, fixed_seam)
    array = zeros(Float64, work_dimension, boundary_dimension, 2, 2)
    same_choices = _edge_choices(same_edges, plan, probabilities)
    cross_choices = _cross_choices(cross_edges, plan, probabilities)
    for left in 1:work_dimension, input in 0:1
        incoming, logical, seam, source =
            _decode_row_work_state(left, cross_count, fixed_seam)
        incoming_parity = _incoming_toggle(incoming, incoming_types)
        for (same_parity, same_logical, same_weight) in same_choices
            output_bit = xor(input, xor(incoming_parity, same_parity))
            partial_source = xor(source, same_parity)
            partial_logical = xor(logical, same_logical)
            for (outgoing, cross_parity, cross_logical, cross_weight) in
                    cross_choices
                xor(partial_source, cross_parity) == source_syndrome ||
                    continue
                right = _row_boundary_state_index(outgoing,
                    xor(partial_logical, cross_logical), seam, cross_count,
                    fixed_seam)
                _add_tensor_entry!(array,
                    (left, right, output_bit + 1, input + 1),
                    same_weight * cross_weight)
            end
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _source_final_tensor(site, left_link, right_link, cross_edges,
        source_syndrome, plan, probabilities, cross_count;
        initial::Bool=false, fixed_seam::Union{Nothing,Int}=nothing)
    output = prime(site)
    boundary_dimension = _row_boundary_state_count(cross_count, fixed_seam)
    cross_choices = _cross_choices(cross_edges, plan, probabilities)
    if initial
        array = zeros(Float64, boundary_dimension, 2, 2)
        seams = isnothing(fixed_seam) ?
            (0:(1 << cross_count)-1) : (fixed_seam:fixed_seam)
        for seam in seams, input in 0:1
            for (outgoing, parity, logical, weight) in cross_choices
                xor(input, parity) == source_syndrome || continue
                right = _row_boundary_state_index(
                    outgoing, logical, seam, cross_count, fixed_seam)
                _add_tensor_entry!(array,
                    (right, 1, input + 1), weight)
            end
        end
        return ITensor(array, right_link, output, site)
    end
    array = zeros(Float64, boundary_dimension, boundary_dimension, 2, 2)
    for left in 1:boundary_dimension, input in 0:1
        incoming, old_logical, seam =
            _decode_row_boundary_state(left, cross_count, fixed_seam)
        for (outgoing, parity, logical, weight) in cross_choices
            xor(input, parity) == source_syndrome || continue
            right = _row_boundary_state_index(outgoing,
                xor(old_logical, logical), seam, cross_count, fixed_seam)
            _add_tensor_entry!(array,
                (left, right, 1, input + 1), weight)
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _logical_closure_tensor(site, left_link, cross_count;
        fixed_seam::Union{Nothing,Int}=nothing)
    boundary_dimension = _row_boundary_state_count(cross_count, fixed_seam)
    array = zeros(Float64, boundary_dimension, 4, 4)
    for left in 1:boundary_dimension, input in 0:3
        incoming, logical, seam = _decode_row_boundary_state(
            left, cross_count, fixed_seam)
        incoming == seam || continue
        array[left, xor(input, logical) + 1, input + 1] = 1.0
    end
    return ITensor(array, left_link, prime(site), site)
end

_fourier_boundary_state_count(cross_count) = 1 << cross_count
_fourier_work_state_count(cross_count) = 2 << cross_count
_fourier_boundary_state_index(incoming) = incoming + 1
_fourier_work_state_index(incoming, source) = 1 + source + 2 * incoming

function _decode_fourier_work_state(index)
    value = index - 1
    return value ÷ 2, value % 2
end

_fourier_character(sector, logical_mask) =
    isodd(count_ones(Int(sector) & Int(logical_mask))) ? -1.0 : 1.0

function _fourier_edge_choices(edge_ids, plan, probabilities, sector)
    choices = Tuple{Int,Float64}[(0, 1.0)]
    for edge_id in edge_ids
        edge = plan.edges[edge_id]
        probability = probabilities[edge_id]
        signed_probability = probability *
            _fourier_character(sector, edge.logical_mask)
        updated = Tuple{Int,Float64}[]
        for (parity, weight) in choices
            iszero((1 - probability) * weight) ||
                push!(updated, (parity, (1 - probability) * weight))
            iszero(signed_probability * weight) ||
                push!(updated,
                    (xor(parity, 1), signed_probability * weight))
        end
        choices = updated
    end
    return choices
end

function _fourier_cross_choices(edge_ids, plan, probabilities, sector)
    choices = Tuple{Int,Int,Float64}[(0, 0, 1.0)]
    for (bit, edge_id) in enumerate(edge_ids)
        edge = plan.edges[edge_id]
        probability = probabilities[edge_id]
        signed_probability = probability *
            _fourier_character(sector, edge.logical_mask)
        updated = Tuple{Int,Int,Float64}[]
        for (mask, parity, weight) in choices
            iszero((1 - probability) * weight) || push!(updated,
                (mask, parity, (1 - probability) * weight))
            iszero(signed_probability * weight) || push!(updated,
                (mask | (1 << (bit - 1)), xor(parity, 1),
                 signed_probability * weight))
        end
        choices = updated
    end
    return choices
end

function _fourier_source_tensor(site, left_link, right_link, cross_count,
        seam; initial::Bool=false)
    output = prime(site)
    work_dimension = _fourier_work_state_count(cross_count)
    if initial
        array = zeros(Float64, work_dimension, 2, 2)
        for input in 0:1
            right = _fourier_work_state_index(seam, input)
            array[right, 1, input + 1] = 1.0
        end
        return ITensor(array, right_link, output, site)
    end
    boundary_dimension = _fourier_boundary_state_count(cross_count)
    array = zeros(Float64, boundary_dimension, work_dimension, 2, 2)
    for incoming in 0:boundary_dimension-1, input in 0:1
        right = _fourier_work_state_index(incoming, input)
        array[incoming + 1, right, 1, input + 1] = 1.0
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _fourier_working_tensor(site, left_link, right_link,
        same_edges, incoming_types, plan, probabilities, cross_count,
        sector)
    output = prime(site)
    work_dimension = _fourier_work_state_count(cross_count)
    array = zeros(Float64, work_dimension, work_dimension, 2, 2)
    choices = _fourier_edge_choices(
        same_edges, plan, probabilities, sector)
    for left in 1:work_dimension, input in 0:1
        incoming, source = _decode_fourier_work_state(left)
        incoming_parity = _incoming_toggle(incoming, incoming_types)
        for (same_parity, weight) in choices
            output_bit = xor(input, xor(incoming_parity, same_parity))
            right = _fourier_work_state_index(
                incoming, xor(source, same_parity))
            _add_tensor_entry!(array,
                (left, right, output_bit + 1, input + 1), weight)
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _fourier_final_working_tensor(site, left_link, right_link,
        same_edges, incoming_types, cross_edges, source_syndrome,
        plan, probabilities, cross_count, sector)
    output = prime(site)
    work_dimension = _fourier_work_state_count(cross_count)
    boundary_dimension = _fourier_boundary_state_count(cross_count)
    array = zeros(Float64, work_dimension, boundary_dimension, 2, 2)
    same_choices = _fourier_edge_choices(
        same_edges, plan, probabilities, sector)
    cross_choices = _fourier_cross_choices(
        cross_edges, plan, probabilities, sector)
    for left in 1:work_dimension, input in 0:1
        incoming, source = _decode_fourier_work_state(left)
        incoming_parity = _incoming_toggle(incoming, incoming_types)
        for (same_parity, same_weight) in same_choices
            output_bit = xor(input, xor(incoming_parity, same_parity))
            partial_source = xor(source, same_parity)
            for (outgoing, cross_parity, cross_weight) in cross_choices
                xor(partial_source, cross_parity) == source_syndrome ||
                    continue
                _add_tensor_entry!(array,
                    (left, outgoing + 1, output_bit + 1, input + 1),
                    same_weight * cross_weight)
            end
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _fourier_source_final_tensor(site, left_link, right_link,
        cross_edges, source_syndrome, plan, probabilities, cross_count,
        sector, seam; initial::Bool=false)
    output = prime(site)
    boundary_dimension = _fourier_boundary_state_count(cross_count)
    choices = _fourier_cross_choices(
        cross_edges, plan, probabilities, sector)
    if initial
        array = zeros(Float64, boundary_dimension, 2, 2)
        for input in 0:1
            for (outgoing, parity, weight) in choices
                xor(input, parity) == source_syndrome || continue
                _add_tensor_entry!(array,
                    (outgoing + 1, 1, input + 1), weight)
            end
        end
        return ITensor(array, right_link, output, site)
    end
    array = zeros(Float64, boundary_dimension, boundary_dimension, 2, 2)
    for incoming in 0:boundary_dimension-1, input in 0:1
        for (outgoing, parity, weight) in choices
            xor(input, parity) == source_syndrome || continue
            _add_tensor_entry!(array,
                (incoming + 1, outgoing + 1, 1, input + 1), weight)
        end
    end
    return ITensor(array, left_link, right_link, output, site)
end

function _fourier_closure_tensor(site, left_link, cross_count, seam)
    boundary_dimension = _fourier_boundary_state_count(cross_count)
    array = zeros(Float64, boundary_dimension, 1, 1)
    array[seam + 1, 1, 1] = 1.0
    return ITensor(array, left_link, prime(site), site)
end

"""
    compile_row_mpo(decoder, strip_id, syndrome)

Compile one conditioned row into an open-chain MPO.  Periodicity is explicit:
the first column sums a seam assignment, the automaton carries it across the
row, and the logical closure tensor requires the last outgoing cross bits to
equal that assignment.
"""
function compile_row_mpo(decoder::RowStripMPSDecoder, strip_id::Integer,
        syndrome; seam::Union{Nothing,Integer}=nothing)
    plan = decoder.plan
    1 <= strip_id <= length(plan.strips) ||
        throw(BoundsError(plan.strips, strip_id))
    observed = _validate_syndrome(plan, syndrome; require_even=false)
    strip = plan.strips[Int(strip_id)]
    sources, edge_tables, same_offsets, cross_offsets =
        _edge_tables(plan, strip)
    cross_count = length(cross_offsets)
    cross_count <= 2 || error(
        "row MPO automaton supports at most two forward cross-column edges")
    fixed_seam = isnothing(seam) ? nothing : Int(seam)
    isnothing(fixed_seam) ||
        (0 <= fixed_seam < (1 << cross_count)) || throw(ArgumentError(
            "seam must be in 0:$(1 << cross_count)-1"))
    boundary_dimension = _row_boundary_state_count(cross_count, fixed_seam)
    work_dimension = _row_work_state_count(cross_count, fixed_seam)

    before = [slot.detector for slot in strip.physical_before]
    after = [slot.detector for slot in strip.physical_after]
    union_detectors = _physical_detector_order(
        plan, unique(vcat(before, after, strip.detectors)))
    columns = [[detector for detector in union_detectors if
        plan.detector_coordinates[detector].column == column]
        for column in 0:plan.L-1]
    all(!isempty, columns) || error("row MPO column has no detector sites")
    for column in 0:plan.L-1
        first(columns[column + 1]) == sources[column] || error(
            "source detector must be the first active role in its column")
    end

    boundary_links = [Index(boundary_dimension,
        "RowBoundary,strip=$strip_id,column=$column")
        for column in 0:plan.L-1]
    tensors = ITensor[]
    for column in 0:plan.L-1
        detectors = columns[column + 1]
        source = sources[column]
        table = edge_tables[column + 1]
        previous_table = edge_tables[mod(column - 1, plan.L) + 1]
        same_by_target = _same_edges_by_target(plan, table, same_offsets)
        incoming_by_target = _incoming_types_by_target(
            plan, previous_table, cross_offsets)
        isempty(get(incoming_by_target, source, Int[])) || error(
            "a cross-column edge unexpectedly targets the source detector")
        cross_edges = [table[offset] for offset in cross_offsets]
        internal_links = [Index(work_dimension,
            "RowWork,strip=$strip_id,column=$column,position=$position")
            for position in 1:max(0, length(detectors) - 1)]
        left_boundary = column == 0 ? nothing : boundary_links[column]
        right_boundary = boundary_links[column + 1]

        for (position, detector) in enumerate(detectors)
            site = decoder.detector_sites[detector]
            first_site = position == 1
            last_site = position == length(detectors)
            if first_site && last_site
                tensor = _source_final_tensor(site, left_boundary,
                    right_boundary, cross_edges, observed[source], plan,
                    decoder.probabilities, cross_count;
                    initial=column == 0, fixed_seam)
            elseif first_site
                tensor = _source_tensor(site, left_boundary,
                    internal_links[1], cross_count; initial=column == 0,
                    fixed_seam)
            elseif last_site
                tensor = _final_working_tensor(site,
                    internal_links[position - 1], right_boundary, detector,
                    get(same_by_target, detector, Int[]),
                    get(incoming_by_target, detector, Int[]), cross_edges,
                    observed[source], plan, decoder.probabilities,
                    cross_count; fixed_seam)
            else
                tensor = _working_tensor(site,
                    internal_links[position - 1], internal_links[position],
                    detector, get(same_by_target, detector, Int[]),
                    get(incoming_by_target, detector, Int[]), plan,
                    decoder.probabilities, cross_count; fixed_seam)
            end
            push!(tensors, tensor)
        end
    end
    push!(tensors, _logical_closure_tensor(
        decoder.logical_site, boundary_links[end], cross_count; fixed_seam))
    length(tensors) == length(union_detectors) + 1 ||
        error("row MPO tensor count is inconsistent")
    return RowMPO(Int(strip_id), union_detectors,
        ITensorMPS.MPO(tensors), cross_offsets,
        max(boundary_dimension, work_dimension), fixed_seam)
end

"""
    compile_fixed_seam_row_mpos(decoder, strip_id, syndrome)

Compile one row as a sum over fixed periodic-seam assignments.  Fixing the
two seam bits removes them from the finite-state automaton, reducing its
maximum dimension from 128 to 32 for the repaired transfer stencil.
"""
function compile_fixed_seam_row_mpos(decoder::RowStripMPSDecoder,
        strip_id::Integer, syndrome)
    # The repaired transfer has two cross-column edge types in all bulk rows.
    # Derive the count instead of assuming it so late rows remain valid.
    1 <= strip_id <= length(decoder.plan.strips) ||
        throw(BoundsError(decoder.plan.strips, strip_id))
    strip = decoder.plan.strips[Int(strip_id)]
    _, _, _, cross_offsets = _edge_tables(decoder.plan, strip)
    return [compile_row_mpo(decoder, strip_id, syndrome; seam)
            for seam in 0:(1 << length(cross_offsets))-1]
end

function _compile_fourier_row_mpo(decoder::RowStripMPSDecoder,
        strip_id::Integer, syndrome, sector::Integer, seam::Integer)
    plan = decoder.plan
    1 <= strip_id <= length(plan.strips) ||
        throw(BoundsError(plan.strips, strip_id))
    0 <= sector <= 3 || throw(ArgumentError("Fourier sector must be in 0:3"))
    observed = _validate_syndrome(plan, syndrome; require_even=false)
    strip = plan.strips[Int(strip_id)]
    sources, edge_tables, same_offsets, cross_offsets =
        _edge_tables(plan, strip)
    cross_count = length(cross_offsets)
    cross_count <= 2 || error(
        "Fourier row MPO supports at most two forward cross-column edges")
    0 <= seam < (1 << cross_count) || throw(ArgumentError(
        "seam must be in 0:$(1 << cross_count)-1"))
    boundary_dimension = _fourier_boundary_state_count(cross_count)
    work_dimension = _fourier_work_state_count(cross_count)

    before = [slot.detector for slot in strip.physical_before]
    after = [slot.detector for slot in strip.physical_after]
    union_detectors = _physical_detector_order(
        plan, unique(vcat(before, after, strip.detectors)))
    columns = [[detector for detector in union_detectors if
        plan.detector_coordinates[detector].column == column]
        for column in 0:plan.L-1]
    all(!isempty, columns) || error(
        "Fourier row MPO column has no detector sites")
    for column in 0:plan.L-1
        first(columns[column + 1]) == sources[column] || error(
            "source detector must be the first active role in its column")
    end

    boundary_links = [Index(boundary_dimension,
        "FourierRowBoundary,sector=$sector,seam=$seam," *
        "strip=$strip_id,column=$column")
        for column in 0:plan.L-1]
    tensors = ITensor[]
    for column in 0:plan.L-1
        detectors = columns[column + 1]
        source = sources[column]
        table = edge_tables[column + 1]
        previous_table = edge_tables[mod(column - 1, plan.L) + 1]
        same_by_target = _same_edges_by_target(plan, table, same_offsets)
        incoming_by_target = _incoming_types_by_target(
            plan, previous_table, cross_offsets)
        isempty(get(incoming_by_target, source, Int[])) || error(
            "a cross-column edge unexpectedly targets the source detector")
        cross_edges = [table[offset] for offset in cross_offsets]
        internal_links = [Index(work_dimension,
            "FourierRowWork,sector=$sector,seam=$seam," *
            "strip=$strip_id,column=$column,position=$position")
            for position in 1:max(0, length(detectors) - 1)]
        left_boundary = column == 0 ? nothing : boundary_links[column]
        right_boundary = boundary_links[column + 1]

        for (position, detector) in enumerate(detectors)
            site = decoder.detector_sites[detector]
            first_site = position == 1
            last_site = position == length(detectors)
            if first_site && last_site
                column == 0 && cross_count != 0 && error(
                    "a one-site first column cannot inject a nonempty seam")
                tensor = _fourier_source_final_tensor(site, left_boundary,
                    right_boundary, cross_edges, observed[source], plan,
                    decoder.probabilities, cross_count, Int(sector),
                    Int(seam); initial=column == 0)
            elseif first_site
                tensor = _fourier_source_tensor(site, left_boundary,
                    internal_links[1], cross_count, Int(seam);
                    initial=column == 0)
            elseif last_site
                tensor = _fourier_final_working_tensor(site,
                    internal_links[position - 1], right_boundary,
                    get(same_by_target, detector, Int[]),
                    get(incoming_by_target, detector, Int[]), cross_edges,
                    observed[source], plan, decoder.probabilities,
                    cross_count, Int(sector))
            else
                tensor = _fourier_working_tensor(site,
                    internal_links[position - 1], internal_links[position],
                    get(same_by_target, detector, Int[]),
                    get(incoming_by_target, detector, Int[]), plan,
                    decoder.probabilities, cross_count, Int(sector))
            end
            push!(tensors, tensor)
        end
    end
    push!(tensors, _fourier_closure_tensor(decoder.scalar_site,
        boundary_links[end], cross_count, Int(seam)))
    length(tensors) == length(union_detectors) + 1 ||
        error("Fourier row MPO tensor count is inconsistent")
    return FourierRowMPO(Int(strip_id), union_detectors,
        ITensorMPS.MPO(tensors), cross_offsets,
        max(boundary_dimension, work_dimension), Int(seam), Int(sector))
end

"""
    compile_fourier_fixed_seam_row_mpos(decoder, strip_id, syndrome, sector)

Compile the fixed-seam MPO branches for one logical Fourier sector.  The
logical class is represented by the edge character
`(-1)^dot(sector, logical_mask)`, so no logical carrier is present in the
MPS and the largest repaired-stencil automaton bond is eight.
"""
function compile_fourier_fixed_seam_row_mpos(
        decoder::RowStripMPSDecoder, strip_id::Integer, syndrome,
        sector::Integer)
    1 <= strip_id <= length(decoder.plan.strips) ||
        throw(BoundsError(decoder.plan.strips, strip_id))
    0 <= sector <= 3 || throw(ArgumentError("Fourier sector must be in 0:3"))
    strip = decoder.plan.strips[Int(strip_id)]
    _, _, _, cross_offsets = _edge_tables(decoder.plan, strip)
    return [_compile_fourier_row_mpo(decoder, strip_id, syndrome,
                sector, seam)
            for seam in 0:(1 << length(cross_offsets))-1]
end

function _site_order(decoder, detectors)
    vcat(Index[decoder.detector_sites[detector] for detector in detectors],
        Index[decoder.logical_site])
end

function _sector_site_order(decoder, detectors)
    vcat(Index[decoder.detector_sites[detector] for detector in detectors],
        Index[decoder.scalar_site])
end

function _boundary_from_vector(decoder, detectors, vector)
    sites = _site_order(decoder, detectors)
    dimensions = Tuple(vcat(fill(2, length(detectors)), 4))
    length(vector) == prod(dimensions) || throw(DimensionMismatch(
        "boundary vector has the wrong dimension"))
    array = reshape(Float64.(vector), dimensions)
    state = ITensorMPS.MPS(array, sites;
        maxdim=decoder.maxdim, cutoff=decoder.cutoff)
    return _BoundaryMPS(state, collect(detectors))
end

function _boundary_vector(decoder, boundary::_BoundaryMPS)
    sites = _site_order(decoder, boundary.detectors)
    tensor = boundary.state[1]
    for position in 2:length(boundary.state)
        tensor *= boundary.state[position]
    end
    return vec(Array(tensor, sites...))
end

function _boundary_mass(decoder, boundary::_BoundaryMPS)
    value = ITensor(1.0)
    sites = _site_order(decoder, boundary.detectors)
    for (position, site) in enumerate(sites)
        value *= boundary.state[position] *
            ITensor(ones(Float64, dim(site)), site)
    end
    return Float64(value[])
end

function _sector_boundary_from_vector(decoder, detectors, vector)
    sites = _sector_site_order(decoder, detectors)
    dimensions = Tuple(vcat(fill(2, length(detectors)), 1))
    length(vector) == prod(dimensions) || throw(DimensionMismatch(
        "sector boundary vector has the wrong dimension"))
    array = reshape(Float64.(vector), dimensions)
    state = ITensorMPS.MPS(array, sites;
        maxdim=decoder.maxdim, cutoff=decoder.cutoff)
    return _SectorBoundaryMPS(state, collect(detectors))
end

function _sector_boundary_vector(decoder, boundary::_SectorBoundaryMPS)
    sites = _sector_site_order(decoder, boundary.detectors)
    tensor = boundary.state[1]
    for position in 2:length(boundary.state)
        tensor *= boundary.state[position]
    end
    return vec(Array(tensor, sites...))
end

function _boundary_mass(decoder, boundary::_SectorBoundaryMPS)
    value = ITensor(1.0)
    sites = _sector_site_order(decoder, boundary.detectors)
    for (position, site) in enumerate(sites)
        value *= boundary.state[position] *
            ITensor(ones(Float64, dim(site)), site)
    end
    return Float64(value[])
end

function _initial_boundary(decoder)
    _boundary_from_vector(decoder, Int[], [1.0, 0.0, 0.0, 0.0])
end

_initial_sector_boundary(decoder) =
    _sector_boundary_from_vector(decoder, Int[], [1.0])

function _maximum_bond_dimension(state)
    maximum_bond = 1
    for position in 1:length(state)-1
        link = ITensorMPS.linkind(state, position)
        isnothing(link) ||
            (maximum_bond = max(maximum_bond, ITensors.dim(link)))
    end
    return maximum_bond
end

function _insert_zero_site(state, position, site, tag)
    1 <= position <= length(state) ||
        throw(BoundsError(1:length(state), position))
    tensors = ITensor[state[index] for index in eachindex(state)]
    if position == 1
        link = Index(1, "RowPad,$tag")
        tensors[1] *= onehot(link => 1)
        insert!(tensors, 1, onehot(site => 1, link => 1))
    else
        old_link = ITensorMPS.linkind(state, position - 1)
        isnothing(old_link) && error("missing MPS link before inserted site")
        new_link = Index(dim(old_link), "RowPad,$tag")
        tensors[position] = replaceind(
            tensors[position], old_link, new_link)
        insert!(tensors, position,
            delta(old_link, new_link) * onehot(site => 1))
    end
    return ITensorMPS.MPS(tensors)
end

_replace_boundary(::_BoundaryMPS, state, detectors) =
    _BoundaryMPS(state, detectors)
_replace_boundary(::_SectorBoundaryMPS, state, detectors) =
    _SectorBoundaryMPS(state, detectors)

function _pad_boundary(decoder, boundary, union_detectors, strip_id)
    state = copy(boundary.state)
    detectors = copy(boundary.detectors)
    for (position, detector) in enumerate(union_detectors)
        if position <= length(detectors) && detectors[position] == detector
            continue
        end
        detector in detectors && error(
            "row union changed the relative detector order")
        state = _insert_zero_site(state, position,
            decoder.detector_sites[detector],
            "strip=$strip_id,detector=$detector")
        insert!(detectors, position, detector)
    end
    detectors == union_detectors || error(
        "padded detector order differs from row MPO")
    return _replace_boundary(boundary, state, detectors)
end

function _remove_zero_outputs(decoder, boundary, output_detectors)
    state = boundary.state
    detectors = copy(boundary.detectors)
    output_set = Set(output_detectors)
    for position in reverse(eachindex(detectors))
        detectors[position] in output_set && continue
        site = decoder.detector_sites[detectors[position]]
        projected = state[position] * onehot(site => 1)
        position < length(state) || error(
            "logical site must remain to the right of detector sites")
        tensors = ITensor[state[index] for index in eachindex(state)]
        tensors[position + 1] = projected * tensors[position + 1]
        deleteat!(tensors, position)
        state = ITensorMPS.MPS(tensors)
        deleteat!(detectors, position)
    end
    detectors == output_detectors || error(
        "stripped detector order differs from compiled active boundary")
    return _replace_boundary(boundary, state, detectors)
end

function _apply_row_mpo(decoder, boundary, row::RowMPO)
    strip = decoder.plan.strips[row.strip_id]
    padded = _pad_boundary(
        decoder, boundary, row.union_detectors, row.strip_id)
    applied = _apply_row_operator(decoder, padded.state, row)
    output_detectors = [slot.detector for slot in strip.physical_after]
    return _remove_zero_outputs(decoder,
        _BoundaryMPS(applied, row.union_detectors), output_detectors)
end

function _apply_row_operator(decoder, padded_state,
        row::Union{RowMPO,FourierRowMPO})
    first_maxdim = decoder.maxdim > typemax(Int) ÷ 2 ?
        typemax(Int) : 2 * decoder.maxdim
    return ITensorMPS.apply(row.operator, padded_state;
        alg="zipup", maxdim=first_maxdim, cutoff=decoder.cutoff / 10,
        truncate_kwargs=(; maxdim=decoder.maxdim, cutoff=decoder.cutoff))
end

function _apply_fourier_fixed_seam_rows(decoder, boundary,
        rows::Vector{FourierRowMPO})
    isempty(rows) && throw(ArgumentError("at least one seam row is required"))
    strip_id = first(rows).strip_id
    sector = first(rows).sector
    all(row -> row.strip_id == strip_id && row.sector == sector, rows) ||
        error("Fourier row branches belong to different strips or sectors")
    union_detectors = first(rows).union_detectors
    all(row -> row.union_detectors == union_detectors, rows) || error(
        "Fourier row branches use different physical sites")
    padded = _pad_boundary(decoder, boundary, union_detectors, strip_id)
    branches = [_apply_row_operator(decoder, padded.state, row)
                for row in rows]
    merged = ITensorMPS.add(branches...; alg="densitymatrix",
        maxdim=decoder.maxdim, cutoff=decoder.cutoff)
    strip = decoder.plan.strips[strip_id]
    output_detectors = [slot.detector for slot in strip.physical_after]
    return _remove_zero_outputs(decoder,
        _SectorBoundaryMPS(merged, union_detectors), output_detectors)
end

function _apply_fixed_seam_rows(decoder, boundary, rows::Vector{RowMPO})
    isempty(rows) && throw(ArgumentError("at least one seam row is required"))
    strip_id = first(rows).strip_id
    all(row -> row.strip_id == strip_id, rows) || error(
        "fixed-seam row branches belong to different strips")
    union_detectors = first(rows).union_detectors
    all(row -> row.union_detectors == union_detectors, rows) || error(
        "fixed-seam row branches use different physical sites")
    padded = _pad_boundary(decoder, boundary, union_detectors, strip_id)
    branches = [_apply_row_operator(decoder, padded.state, row)
                for row in rows]
    merged = ITensorMPS.add(branches...; alg="densitymatrix",
        maxdim=decoder.maxdim, cutoff=decoder.cutoff)
    strip = decoder.plan.strips[strip_id]
    output_detectors = [slot.detector for slot in strip.physical_after]
    return _remove_zero_outputs(decoder,
        _BoundaryMPS(merged, union_detectors), output_detectors)
end

"""Apply one compiled row MPO to a dense boundary vector (test/reference API)."""
function row_mpo_action(decoder::RowStripMPSDecoder, row::RowMPO, input)
    strip = decoder.plan.strips[row.strip_id]
    input_detectors = [slot.detector for slot in strip.physical_before]
    boundary = _boundary_from_vector(decoder, input_detectors, input)
    output = _apply_row_mpo(decoder, boundary, row)
    return _boundary_vector(decoder, output)
end

function _contract_fourier_sector(decoder::RowStripMPSDecoder, observed,
        sector::Int)
    boundary = _initial_sector_boundary(decoder)
    logscale = 0.0
    maximum_bond = 1
    maximum_automaton_dimension = 1
    mpo_applications = 0
    svd_events = 0
    for strip in decoder.plan.strips
        boundary.detectors ==
            [slot.detector for slot in strip.physical_before] ||
            error("Fourier row boundary order differs from compiled plan")
        rows = compile_fourier_fixed_seam_row_mpos(
            decoder, strip.id, observed, sector)
        boundary = _apply_fourier_fixed_seam_rows(decoder, boundary, rows)
        mpo_applications += length(rows)
        maximum_automaton_dimension = max(
            maximum_automaton_dimension,
            maximum(row.maximum_automaton_dimension for row in rows))
        svd_events += (2 * length(rows) + 1) *
            max(0, length(first(rows).operator) - 1)
        magnitude = ITensorMPS.norm(boundary.state)
        isfinite(magnitude) && magnitude >= 0 || throw(DomainError(
            magnitude, "invalid Fourier-sector MPS norm at " *
                "sector $sector strip $(strip.id)"))
        if iszero(magnitude)
            return (sign=0, logabs=-Inf, maximum_bond,
                maximum_automaton_dimension, mpo_applications, svd_events)
        end
        boundary.state[1] /= magnitude
        logscale += log(magnitude)
        maximum_bond = max(maximum_bond,
            _maximum_bond_dimension(boundary.state))
    end

    isempty(boundary.detectors) || error(
        "final Fourier-sector boundary is not empty")
    scalar = only(_sector_boundary_vector(decoder, boundary))
    isfinite(scalar) || throw(DomainError(scalar,
        "nonfinite final Fourier-sector scalar"))
    return (sign=iszero(scalar) ? 0 : Int(sign(scalar)),
        logabs=iszero(scalar) ? -Inf : logscale + log(abs(scalar)),
        maximum_bond, maximum_automaton_dimension, mpo_applications,
        svd_events)
end

function _decode_fourier_sectors(decoder::RowStripMPSDecoder, observed;
        sector_parallel::Bool=false)
    results = Vector{NamedTuple}(undef, 4)
    parallel_used = sector_parallel && Threads.nthreads() > 1
    if parallel_used
        @sync for sector in 0:3
            Threads.@spawn results[sector + 1] =
                _contract_fourier_sector(decoder, observed, sector)
        end
    else
        for sector in 0:3
            results[sector + 1] =
                _contract_fourier_sector(decoder, observed, sector)
        end
    end
    sector_signs = [result.sign for result in results]
    sector_logabs = [result.logabs for result in results]
    maximum_bond = maximum(result.maximum_bond for result in results)
    maximum_automaton_dimension = maximum(
        result.maximum_automaton_dimension for result in results)
    mpo_applications = sum(result.mpo_applications for result in results)
    svd_events = sum(result.svd_events for result in results)

    sector_signs[1] == 1 && isfinite(sector_logabs[1]) ||
        throw(DomainError((sector_signs[1], sector_logabs[1]),
            "the positive Fourier sector must have positive finite mass"))
    ratios = [sector_signs[index] == 0 ? 0.0 :
        sector_signs[index] * exp(sector_logabs[index] - sector_logabs[1])
        for index in 1:4]
    raw_weights = [0.25 * sum(
        _fourier_character(sector, logical) * ratios[sector + 1]
        for sector in 0:3) for logical in 0:3]
    negative_tolerance = 1e-10
    ratio_violation = max(0.0, maximum(abs, ratios) - 1.0)
    roundoff_correction_l1 = sum(
        (-weight for weight in raw_weights
         if -negative_tolerance <= weight < 0); init=0.0)
    valid = all(isfinite, raw_weights) &&
        minimum(raw_weights) >= -negative_tolerance &&
        ratio_violation <= negative_tolerance
    weights = [weight < 0 && weight >= -negative_tolerance ? 0.0 : weight
               for weight in raw_weights]
    total = sum(weights)
    valid &= isfinite(total) && total > 0
    probabilities = valid ? weights ./ total : fill(NaN, 4)
    diagnostics = (
        algorithm="fourier_sector_fixed_seam_row_mpo_boundary_mps",
        sector_count=4,
        sector_parallel=parallel_used,
        julia_threads=Threads.nthreads(),
        sector_signs,
        sector_logabs,
        row_updates=4 * length(decoder.plan.strips),
        svd_events,
        transition_count=0,
        mpo_applications,
        maximum_automaton_dimension,
        actual_maxbond=maximum_bond,
        maxdim=decoder.maxdim,
        cutoff=decoder.cutoff,
        exact_bond_guard=decoder.exact_bond_guard,
        # The final active-boundary rank does not bound the transient
        # MPO×MPS rank encountered by zip-up. Without explicit discarded-rank
        # diagnostics this backend cannot certify an untruncated contraction.
        untruncated=false,
        boundary_rank_guard_satisfied=decoder.cutoff == 0.0 &&
            decoder.maxdim >= decoder.exact_bond_guard,
        maximum_active_boundary=decoder.plan.maximum_active_boundary,
        roundoff_correction_l1,
        maximum_fourier_ratio=maximum(abs, ratios),
        ratio_violation,
    )
    return (valid, probabilities,
        logical_mask=valid ? argmax(probabilities) - 1 : -1,
        log_syndrome=sector_logabs[1], raw_probabilities=raw_weights,
        diagnostics)
end

function decode_row_strip_mps!(decoder::RowStripMPSDecoder, syndrome;
        algorithm::Symbol=:row_mpo, sector_parallel::Bool=false)
    algorithm in (:dense, :row_mpo, :row_mpo_fixed_seams,
                  :row_mpo_fourier_sectors) ||
        throw(ArgumentError(
            "algorithm must be :dense, :row_mpo, " *
            ":row_mpo_fixed_seams, or :row_mpo_fourier_sectors"))
    observed = _validate_syndrome(decoder.plan, syndrome)
    algorithm === :row_mpo_fourier_sectors &&
        return _decode_fourier_sectors(
            decoder, observed; sector_parallel)
    sector_parallel && throw(ArgumentError(
        "sector_parallel is only valid for :row_mpo_fourier_sectors"))
    boundary = _initial_boundary(decoder)
    log_syndrome = 0.0
    maximum_bond = 1
    svd_events = 0
    transition_count = 0
    mpo_applications = 0
    maximum_automaton_dimension = 1

    for strip in decoder.plan.strips
        boundary.detectors ==
            [slot.detector for slot in strip.physical_before] ||
            error("row boundary order differs from the compiled plan")
        if algorithm === :dense
            row = compile_dense_row_map(decoder.plan, strip.id,
                decoder.probabilities, observed)
            output = row.matrix * _boundary_vector(decoder, boundary)
            boundary = _boundary_from_vector(
                decoder, row.output_detectors, output)
            transition_count += row.transition_count
            svd_events += max(0, length(boundary.state) - 1)
        elseif algorithm === :row_mpo
            row = compile_row_mpo(decoder, strip.id, observed)
            boundary = _apply_row_mpo(decoder, boundary, row)
            mpo_applications += 1
            maximum_automaton_dimension = max(
                maximum_automaton_dimension,
                row.maximum_automaton_dimension)
            svd_events += 2 * max(0, length(row.operator) - 1)
        else
            rows = compile_fixed_seam_row_mpos(
                decoder, strip.id, observed)
            boundary = _apply_fixed_seam_rows(decoder, boundary, rows)
            mpo_applications += length(rows)
            maximum_automaton_dimension = max(
                maximum_automaton_dimension,
                maximum(row.maximum_automaton_dimension for row in rows))
            svd_events += (2 * length(rows) + 1) *
                max(0, length(first(rows).operator) - 1)
        end
        mass = _boundary_mass(decoder, boundary)
        isfinite(mass) && mass > 0 || throw(DomainError(mass,
            "nonpositive or nonfinite conditioned mass at strip $(strip.id)"))
        boundary.state[1] /= mass
        log_syndrome += log(mass)
        maximum_bond = max(maximum_bond,
            _maximum_bond_dimension(boundary.state))
    end

    isempty(boundary.detectors) ||
        error("final row-strip boundary is not empty")
    weights = _boundary_vector(decoder, boundary)
    total = sum(weights)
    isfinite(total) && total > 0 ||
        throw(DomainError(total, "invalid final logical mass"))
    probabilities = weights ./ total
    log_syndrome += log(total)
    valid = all(isfinite, probabilities) && minimum(probabilities) >= -1e-10
    diagnostics = (
        algorithm=algorithm === :dense ?
            "dense_row_map_boundary_mps_semantic_prototype" :
            algorithm === :row_mpo ?
                "open_chain_row_mpo_boundary_mps" :
                "fixed_seam_open_chain_row_mpo_boundary_mps",
        row_updates=length(decoder.plan.strips),
        svd_events,
        transition_count,
        mpo_applications,
        maximum_automaton_dimension,
        actual_maxbond=maximum_bond,
        maxdim=decoder.maxdim,
        cutoff=decoder.cutoff,
        exact_bond_guard=decoder.exact_bond_guard,
        untruncated=algorithm === :dense && decoder.cutoff == 0.0 &&
            decoder.maxdim >= decoder.exact_bond_guard,
        maximum_active_boundary=decoder.plan.maximum_active_boundary,
    )
    return (valid, probabilities,
        logical_mask=valid ? argmax(probabilities) - 1 : -1,
        log_syndrome, raw_probabilities=probabilities, diagnostics)
end

"""
Build and run the experimental row-strip boundary-MPS decoder.

For `algorithm=:row_mpo_fourier_sectors`, `sector_parallel=true` schedules
the four independent Fourier sectors on Julia threads.  Run Julia with at
least four threads and keep BLAS at one thread to avoid oversubscription.
"""
function decode_row_strip_mps(topology, L::Integer, p::Real, syndrome;
        maxdim::Integer=4096, cutoff::Real=0.0,
        algorithm::Symbol=:row_mpo, sector_parallel::Bool=false)
    decoder = build_row_strip_mps_decoder(topology, L, p;
        maxdim, cutoff)
    return decode_row_strip_mps!(decoder, syndrome;
        algorithm, sector_parallel)
end

end # module
