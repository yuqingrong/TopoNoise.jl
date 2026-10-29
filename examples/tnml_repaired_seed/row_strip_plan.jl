module TNMLRowStripPlan

export DetectorCoordinate, RowPlane, DetectorLayout, BoundarySlot, EdgeGeometry,
       RowStrip, RowStripPlan, repaired_seed_transfer_layout,
       compile_repaired_seed_transfer_plan, audit_row_strip_plan

"""One final-star detector in the repaired-seed transfer factor graph."""
struct DetectorCoordinate
    id::Int
    kind::Symbol
    plane::Int
    time::Int
    row::Int
    column::Int
end

Base.:(==)(a::DetectorCoordinate, b::DetectorCoordinate) =
    a.id == b.id && a.kind === b.kind && a.plane == b.plane &&
    a.time == b.time && a.row == b.row && a.column == b.column

"""A plane of detector rows in plane-major sweep order."""
struct RowPlane
    id::Int
    kind::Symbol
    time::Int
    detector_ids::Vector{Int}
    strip_ids::Vector{Int}
end

"""The fixed detector geometry, independent of faults and probabilities."""
struct DetectorLayout
    L::Int
    detectors::Vector{DetectorCoordinate}
    planes::Vector{RowPlane}
end

"""A detector carried across a row cut, with its position at that cut."""
struct BoundarySlot
    slot::Int
    detector::Int
    coordinate::DetectorCoordinate
end

"""Geometry and ownership of one effective independent edge variable."""
struct EdgeGeometry
    id::Int
    endpoints::NTuple{2,Int}
    logical_mask::UInt8
    owner_strip::Int
    offset::NTuple{3,Int}
end

"""
One row update: introduce owned edges, impose syndrome, then retire detectors.

`active_*` is the canonical detector-ID order used for identity and audits.
`physical_*` groups roles by periodic column using key `(column, detector_id)`
for a numerical MPS. `canonical_to_physical_*[i]` gives the physical slot of
canonical slot `i`.
"""
struct RowStrip
    id::Int
    plane::Int
    row::Int
    detectors::Vector{Int}
    owned_edges::Vector{Int}
    retired_detectors::Vector{Int}
    active_before::Vector{BoundarySlot}
    active_after::Vector{BoundarySlot}
    physical_before::Vector{BoundarySlot}
    physical_after::Vector{BoundarySlot}
    canonical_to_physical_before::Vector{Int}
    canonical_to_physical_after::Vector{Int}
end

"""
Exact geometry plan for a plane-major, row-major transfer-network sweep.

After strip `k`, the active slots are exactly the future endpoints of edges
introduced by strips `1:k` whose detector has not yet been retired. The
two-bit logical class is deliberately absent from these detector slots: a
numerical implementation keeps it as one separate four-state open index.
"""
struct RowStripPlan
    L::Int
    detector_coordinates::Vector{DetectorCoordinate}
    planes::Vector{RowPlane}
    edges::Vector{EdgeGeometry}
    strips::Vector{RowStrip}
    detector_to_strip::Vector{Int}
    edge_owner::Vector{Int}
    allowed_offsets::Set{NTuple{3,Int}}
    maximum_active_boundary::Int
end

function _checked_L(L::Integer)
    2 <= L <= typemax(Int) ||
        throw(ArgumentError("L must be at least two and fit Int"))
    (big(L) + 2) * L * L <= typemax(Int) ||
        throw(ArgumentError("detector count must fit Int"))
    return Int(L)
end

function _plane_kind(plane::Int, L::Int)
    plane == 0 && return :seed_white
    plane == 1 && return :initial_transfer
    plane == L + 1 && return :final_frontier
    return :bulk
end

"""
    repaired_seed_transfer_layout(L)

Lay out the `(L+2)L^2` final-star detectors used by the transfer component
after the seed repair. Plane zero contains the initial white seed stars.
Planes `1:L+1` contain circuit times `0:L`, with the final frontier last.
Within every plane, detectors are ordered by row `u`, then periodic column `v`.
"""
function repaired_seed_transfer_layout(L::Integer)
    size = _checked_L(L)
    detectors = DetectorCoordinate[]
    planes = RowPlane[]
    for plane in 0:size+1
        kind = _plane_kind(plane, size)
        time = plane - 1
        first_detector = length(detectors) + 1
        for row in 0:size-1, column in 0:size-1
            push!(detectors, DetectorCoordinate(length(detectors) + 1, kind,
                plane, time, row, column))
        end
        last_detector = length(detectors)
        strip_ids = [plane * size + row + 1 for row in 0:size-1]
        push!(planes, RowPlane(plane + 1, kind, time,
            collect(first_detector:last_detector), strip_ids))
    end
    return DetectorLayout(size, detectors, planes)
end

# The diagonal H+6CNOT circuit only connects these relative final-star
# coordinates. Transverse offsets are forward periodic differences.
function _repaired_transfer_offsets()
    Set{NTuple{3,Int}}([
        (1, 0, 0), (1, 0, 1), (1, 1, 0),
        (2, 1, 0), (2, 1, 1),
    ])
end

function _boundary_slots(detectors, coordinates)
    [BoundarySlot(slot, detector, coordinates[detector])
     for (slot, detector) in enumerate(sort!(collect(detectors)))]
end

# Numerical tensors should group all roles belonging to one periodic column.
# The canonical boundary above stays detector-sorted for identity and audits;
# this separate order prevents detector-ID plane blocks from separating local
# v-column factors by O(L) sites.
function _physical_slots(detectors, coordinates)
    ordered = sort!(collect(detectors); by=detector -> begin
        coordinate = coordinates[detector]
        (coordinate.column, coordinate.id)
    end)
    [BoundarySlot(slot, detector, coordinates[detector])
     for (slot, detector) in enumerate(ordered)]
end

function _canonical_to_physical(canonical, physical)
    positions = Dict(slot.detector => slot.slot for slot in physical)
    [positions[slot.detector] for slot in canonical]
end

"""Compile directly from an object exposing effective endpoints and masks."""
function compile_repaired_seed_transfer_plan(topology, L::Integer)
    hasproperty(topology, :detector_endpoints) &&
        hasproperty(topology, :logical_masks) || throw(ArgumentError(
        "topology must expose detector_endpoints and logical_masks"))
    return compile_repaired_seed_transfer_plan(
        getproperty(topology, :detector_endpoints),
        getproperty(topology, :logical_masks), L)
end

"""
    compile_repaired_seed_transfer_plan(endpoints, logical_masks, L)

Assign each effective edge to the strip containing its earlier endpoint. A
strip introduces all of those edge factors before imposing and retiring its
entire detector row. No probability or observed syndrome enters this plan.
"""
function compile_repaired_seed_transfer_plan(endpoints::AbstractVector,
        logical_masks::AbstractVector, L::Integer)
    layout = repaired_seed_transfer_layout(L)
    size = layout.L
    length(endpoints) == length(logical_masks) ||
        throw(DimensionMismatch("edge endpoints and logical masks differ in length"))
    detector_count = length(layout.detectors)
    strip_count = (size + 2) * size
    detector_to_strip = Vector{Int}(undef, detector_count)
    strip_detectors = [Int[] for _ in 1:strip_count]
    for coordinate in layout.detectors
        strip = coordinate.plane * size + coordinate.row + 1
        detector_to_strip[coordinate.id] = strip
        push!(strip_detectors[strip], coordinate.id)
    end

    allowed_offsets = _repaired_transfer_offsets()
    edge_owner = Vector{Int}(undef, length(endpoints))
    owned_edges = [Int[] for _ in 1:strip_count]
    edges = EdgeGeometry[]
    signatures = Set{Tuple{NTuple{2,Int},UInt8}}()
    for (edge_id, (raw_endpoint, raw_mask)) in
            enumerate(zip(endpoints, logical_masks))
        length(raw_endpoint) == 2 ||
            throw(ArgumentError("every effective edge needs two endpoints"))
        source, destination = Int(raw_endpoint[1]), Int(raw_endpoint[2])
        1 <= source < destination <= detector_count || throw(ArgumentError(
            "effective-edge endpoints must be ordered, distinct, and in bounds"))
        raw_mask isa Integer && 0 <= raw_mask <= 3 ||
            throw(ArgumentError("logical masks must lie in 0:3"))
        signature = ((source, destination), UInt8(raw_mask))
        signature in signatures && throw(ArgumentError(
            "effective endpoint and logical-mask signatures must be unique"))
        push!(signatures, signature)
        first = layout.detectors[source]
        second = layout.detectors[destination]
        offset = (second.plane - first.plane,
                  mod(second.row - first.row, size),
                  mod(second.column - first.column, size))
        offset in allowed_offsets || throw(ArgumentError(
            "edge $edge_id has circuit-forbidden detector offset $offset"))
        owner = detector_to_strip[source]
        edge_owner[edge_id] = owner
        push!(owned_edges[owner], edge_id)
        push!(edges, EdgeGeometry(edge_id, (source, destination),
            UInt8(raw_mask), owner, offset))
    end

    strips = RowStrip[]
    active = Set{Int}()
    for strip_id in 1:strip_count
        before = _boundary_slots(active, layout.detectors)
        physical_before = _physical_slots(active, layout.detectors)
        for edge_id in owned_edges[strip_id]
            push!(active, edges[edge_id].endpoints[2])
        end
        detectors = strip_detectors[strip_id]
        foreach(detector -> delete!(active, detector), detectors)
        after = _boundary_slots(active, layout.detectors)
        physical_after = _physical_slots(active, layout.detectors)
        before_permutation = _canonical_to_physical(before, physical_before)
        after_permutation = _canonical_to_physical(after, physical_after)
        plane, row = divrem(strip_id - 1, size)
        push!(strips, RowStrip(strip_id, plane, row, copy(detectors),
            copy(owned_edges[strip_id]), copy(detectors), before, after,
            physical_before, physical_after, before_permutation,
            after_permutation))
    end
    maximum_active = maximum(length(strip.active_after) for strip in strips;
        init=0)
    plan = RowStripPlan(size, layout.detectors, layout.planes, edges, strips,
        detector_to_strip, edge_owner, allowed_offsets, maximum_active)
    audit_row_strip_plan(plan)
    return plan
end

"""Recompute all ownership, retirement, geometry, and boundary invariants."""
function audit_row_strip_plan(plan::RowStripPlan)
    size = plan.L
    detector_count = length(plan.detector_coordinates)
    edge_count = length(plan.edges)
    strip_count = length(plan.strips)
    detector_count == (size + 2) * size^2 ||
        error("detector layout has the wrong size")
    strip_count == (size + 2) * size ||
        error("row-strip layout has the wrong size")
    length(plan.detector_to_strip) == detector_count ||
        error("detector-to-strip map has the wrong size")
    length(plan.edge_owner) == edge_count ||
        error("edge-owner map has the wrong size")
    plan.allowed_offsets == _repaired_transfer_offsets() ||
        error("allowed circuit offsets have changed")
    expected_layout = repaired_seed_transfer_layout(size)
    plan.detector_coordinates == expected_layout.detectors ||
        error("detector coordinates do not match the circuit layout")
    length(plan.planes) == length(expected_layout.planes) ||
        error("detector plane count is inconsistent")
    for (plane, expected) in zip(plan.planes, expected_layout.planes)
        plane.id == expected.id && plane.kind === expected.kind &&
            plane.time == expected.time &&
            plane.detector_ids == expected.detector_ids &&
            plane.strip_ids == expected.strip_ids ||
            error("detector plane layout is inconsistent")
    end

    retired = zeros(Int, detector_count)
    owned = zeros(Int, edge_count)
    active = Set{Int}()
    for (strip_id, strip) in enumerate(plan.strips)
        strip.id == strip_id || error("strip IDs are not consecutive")
        strip.plane * size + strip.row + 1 == strip_id ||
            error("strip coordinate does not match its ID")
        before = sort!(collect(active))
        [slot.detector for slot in strip.active_before] == before ||
            error("active-before boundary is inconsistent")
        [slot.slot for slot in strip.active_before] ==
            collect(1:length(before)) || error("active-before slots are not dense")
        all(slot.coordinate == plan.detector_coordinates[slot.detector]
            for slot in strip.active_before) ||
            error("active-before slot coordinates are inconsistent")
        physical_before = sort(before; by=detector -> begin
            coordinate = plan.detector_coordinates[detector]
            (coordinate.column, coordinate.id)
        end)
        [slot.detector for slot in strip.physical_before] == physical_before ||
            error("physical-before boundary order is inconsistent")
        [slot.slot for slot in strip.physical_before] ==
            collect(1:length(physical_before)) ||
            error("physical-before slots are not dense")
        all(slot.coordinate == plan.detector_coordinates[slot.detector]
            for slot in strip.physical_before) ||
            error("physical-before slot coordinates are inconsistent")
        strip.canonical_to_physical_before ==
            _canonical_to_physical(strip.active_before, strip.physical_before) ||
            error("active-before physical permutation is inconsistent")
        strip.retired_detectors == strip.detectors ||
            error("a strip must retire exactly its detector row")
        for detector in strip.detectors
            1 <= detector <= detector_count || error("detector is out of bounds")
            plan.detector_to_strip[detector] == strip_id ||
                error("detector belongs to the wrong strip")
            retired[detector] += 1
        end
        for edge_id in strip.owned_edges
            1 <= edge_id <= edge_count || error("edge is out of bounds")
            edge = plan.edges[edge_id]
            edge.id == edge_id || error("edge IDs are not consecutive")
            edge.owner_strip == strip_id == plan.edge_owner[edge_id] ||
                error("edge owner is inconsistent")
            plan.detector_to_strip[edge.endpoints[1]] == strip_id ||
                error("edge is not owned by its source detector row")
            plan.detector_to_strip[edge.endpoints[2]] > strip_id ||
                error("edge destination is not in a future strip")
            edge.offset in plan.allowed_offsets ||
                error("edge has a forbidden geometry offset")
            source = plan.detector_coordinates[edge.endpoints[1]]
            destination = plan.detector_coordinates[edge.endpoints[2]]
            expected_offset = (destination.plane - source.plane,
                mod(destination.row - source.row, size),
                mod(destination.column - source.column, size))
            edge.offset == expected_offset ||
                error("edge offset does not match its detector coordinates")
            owned[edge_id] += 1
            push!(active, edge.endpoints[2])
        end
        foreach(detector -> delete!(active, detector), strip.detectors)
        after = sort!(collect(active))
        [slot.detector for slot in strip.active_after] == after ||
            error("active-after boundary is inconsistent")
        [slot.slot for slot in strip.active_after] ==
            collect(1:length(after)) || error("active-after slots are not dense")
        all(slot.coordinate == plan.detector_coordinates[slot.detector]
            for slot in strip.active_after) ||
            error("active-after slot coordinates are inconsistent")
        physical_after = sort(after; by=detector -> begin
            coordinate = plan.detector_coordinates[detector]
            (coordinate.column, coordinate.id)
        end)
        [slot.detector for slot in strip.physical_after] == physical_after ||
            error("physical-after boundary order is inconsistent")
        [slot.slot for slot in strip.physical_after] ==
            collect(1:length(physical_after)) ||
            error("physical-after slots are not dense")
        all(slot.coordinate == plan.detector_coordinates[slot.detector]
            for slot in strip.physical_after) ||
            error("physical-after slot coordinates are inconsistent")
        strip.canonical_to_physical_after ==
            _canonical_to_physical(strip.active_after, strip.physical_after) ||
            error("active-after physical permutation is inconsistent")
    end
    all(==(1), retired) || error("every detector must retire exactly once")
    all(==(1), owned) || error("every edge must be owned exactly once")
    isempty(active) || error("final active boundary is not empty")
    maximum_active = maximum(length(strip.active_after) for strip in plan.strips;
        init=0)
    maximum_active == plan.maximum_active_boundary ||
        error("recorded maximum active boundary is inconsistent")
    maximum_plane_span = maximum((edge.offset[1] for edge in plan.edges);
        init=0)
    observed_offsets = Set(edge.offset for edge in plan.edges)
    return (valid=true, detector_count, edge_count, strip_count,
            maximum_plane_span, maximum_active_boundary=maximum_active,
            allowed_offsets=copy(plan.allowed_offsets), observed_offsets)
end

end # module
