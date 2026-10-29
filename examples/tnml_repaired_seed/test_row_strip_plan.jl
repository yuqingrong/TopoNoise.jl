using Test, TopoNoise

include("model.jl")
include("row_strip_plan.jl")

using .RepairedSeedTN
using .TNMLRowStripPlan

function direct_factor_distribution(detector_count, endpoints, logical_masks,
        probabilities)
    length(endpoints) == length(logical_masks) == length(probabilities) ||
        throw(DimensionMismatch("edge arrays differ in length"))
    states = zeros(Float64, 1 << (detector_count + 2))
    scratch = similar(states)
    states[1] = 1.0
    for (edge, logical, probability) in
            zip(endpoints, logical_masks, probabilities)
        u, v = edge
        signature = (1 << (u - 1)) | (1 << (v - 1)) |
                    (Int(logical) << detector_count)
        @inbounds for state in 0:length(states)-1
            scratch[state + 1] =
                (1 - probability) * states[state + 1] +
                probability * states[xor(state, signature) + 1]
        end
        states, scratch = scratch, states
    end
    return states
end

function direct_factor_weights(states, detector_count, syndrome)
    syndrome_mask = sum(Int(bit) << (detector - 1)
        for (detector, bit) in enumerate(syndrome); init=0)
    return [states[(syndrome_mask | (logical << detector_count)) + 1]
            for logical in 0:3]
end

function remap_mask(mask::BigInt, old_detectors, new_detectors)
    positions = Dict(detector => position
        for (position, detector) in enumerate(new_detectors))
    result = big(0)
    for (old_position, detector) in enumerate(old_detectors)
        iszero(mask & (big(1) << (old_position - 1))) && continue
        haskey(positions, detector) || continue
        result |= big(1) << (positions[detector] - 1)
    end
    return result
end

"""Independent exact frontier contraction driven only by the compiled plan."""
function row_plan_weights(plan, probabilities, syndrome)
    length(probabilities) == length(plan.edges) ||
        throw(DimensionMismatch("incorrect edge probability count"))
    length(syndrome) == length(plan.detector_coordinates) ||
        throw(DimensionMismatch("incorrect syndrome length"))
    states = Dict((big(0), UInt8(0)) => 1.0)
    active = Int[]
    for strip in plan.strips
        @test active == [slot.detector for slot in strip.active_before]
        workspace = sort!(unique(vcat(active, strip.detectors,
            [endpoint for edge_id in strip.owned_edges
                for endpoint in plan.edges[edge_id].endpoints])))
        expanded = Dict{Tuple{BigInt,UInt8},Float64}()
        for ((mask, logical), weight) in states
            key = (remap_mask(mask, active, workspace), logical)
            expanded[key] = get(expanded, key, 0.0) + weight
        end
        positions = Dict(detector => position
            for (position, detector) in enumerate(workspace))
        for edge_id in strip.owned_edges
            edge = plan.edges[edge_id]
            probability = probabilities[edge_id]
            delta = (big(1) << (positions[edge.endpoints[1]] - 1)) |
                    (big(1) << (positions[edge.endpoints[2]] - 1))
            next = Dict{Tuple{BigInt,UInt8},Float64}()
            for ((mask, logical), weight) in expanded
                absent = (mask, logical)
                present = (xor(mask, delta), xor(logical, edge.logical_mask))
                next[absent] = get(next, absent, 0.0) +
                    (1 - probability) * weight
                next[present] = get(next, present, 0.0) + probability * weight
            end
            expanded = next
        end
        after = [slot.detector for slot in strip.active_after]
        states = Dict{Tuple{BigInt,UInt8},Float64}()
        for ((mask, logical), weight) in expanded
            all(detector ->
                !iszero(mask & (big(1) << (positions[detector] - 1))) ==
                    Bool(syndrome[detector]), strip.retired_detectors) || continue
            key = (remap_mask(mask, workspace, after), logical)
            states[key] = get(states, key, 0.0) + weight
        end
        active = after
    end
    @test isempty(active)
    return [get(states, (big(0), UInt8(logical)), 0.0) for logical in 0:3]
end

@testset "repaired-seed transfer detector layout" begin
    layout = repaired_seed_transfer_layout(2)
    @test length(layout.detectors) == 16
    @test length(layout.planes) == 4
    @test layout.detectors[1] == DetectorCoordinate(
        1, :seed_white, 0, -1, 0, 0)
    @test layout.detectors[4] == DetectorCoordinate(
        4, :seed_white, 0, -1, 1, 1)
    @test layout.detectors[5] == DetectorCoordinate(
        5, :initial_transfer, 1, 0, 0, 0)
    @test layout.detectors[16] == DetectorCoordinate(
        16, :final_frontier, 3, 2, 1, 1)
    @test [plane.kind for plane in layout.planes] ==
        [:seed_white, :initial_transfer, :bulk, :final_frontier]
    @test all(length(plane.strip_ids) == 2 for plane in layout.planes)
end

@testset "row-strip ownership and active-boundary audit" begin
    context = build_context(2)
    topology = context.transfer_topology
    plan = compile_repaired_seed_transfer_plan(topology, 2)
    audit = audit_row_strip_plan(plan)

    @test audit.valid
    @test audit.detector_count == 16
    @test audit.edge_count == 52
    @test audit.strip_count == 8
    @test audit.maximum_plane_span == 2
    @test audit.allowed_offsets == Set([
        (1, 0, 0), (1, 0, 1), (1, 1, 0),
        (2, 1, 0), (2, 1, 1)])
    @test [edge.endpoints for edge in plan.edges] ==
        topology.detector_endpoints
    @test [edge.logical_mask for edge in plan.edges] ==
        topology.logical_masks
    @test sort(vcat((strip.detectors for strip in plan.strips)...)) == 1:16
    @test sort(vcat((strip.retired_detectors for strip in plan.strips)...)) == 1:16
    @test sort(vcat((strip.owned_edges for strip in plan.strips)...)) == 1:52
    @test isempty(plan.strips[1].active_before)
    @test isempty(plan.strips[end].active_after)
    @test plan.maximum_active_boundary == maximum(length(strip.active_after)
        for strip in plan.strips)

    for strip in plan.strips
        @test [slot.slot for slot in strip.active_before] ==
            collect(1:length(strip.active_before))
        @test [slot.slot for slot in strip.active_after] ==
            collect(1:length(strip.active_after))
        @test issorted(slot.detector for slot in strip.active_before)
        @test issorted(slot.detector for slot in strip.active_after)
        for (canonical, physical, permutation) in
                ((strip.active_before, strip.physical_before,
                  strip.canonical_to_physical_before),
                 (strip.active_after, strip.physical_after,
                  strip.canonical_to_physical_after))
            @test sort([slot.detector for slot in physical]) ==
                [slot.detector for slot in canonical]
            @test [slot.slot for slot in physical] == collect(1:length(physical))
            @test [slot.detector for slot in physical] == sort(
                [slot.detector for slot in canonical]; by=detector -> begin
                    coordinate = plan.detector_coordinates[detector]
                    (coordinate.column, coordinate.id)
                end)
            @test permutation == [findfirst(slot ->
                slot.detector == canonical_slot.detector, physical)
                for canonical_slot in canonical]
        end
        for edge_id in strip.owned_edges
            edge = plan.edges[edge_id]
            @test edge.owner_strip == strip.id
            @test plan.detector_to_strip[edge.endpoints[1]] == strip.id
            @test edge.offset in audit.allowed_offsets
        end
    end
end

@testset "larger repaired-transfer graph keeps the local stencil" begin
    context = build_context(3)
    plan = compile_repaired_seed_transfer_plan(context.transfer_topology, 3)
    audit = audit_row_strip_plan(plan)
    @test audit.detector_count == 45
    @test audit.edge_count == 162
    @test audit.strip_count == 15
    @test audit.maximum_plane_span == 2
    @test audit.observed_offsets == audit.allowed_offsets
    @test isempty(plan.strips[end].active_after)
end

@testset "row plan preserves exact four-logical factor weights" begin
    context = build_context(2)
    topology = context.transfer_topology
    plan = compile_repaired_seed_transfer_plan(topology, 2)
    for p in (0.02, 0.1)
        probabilities = toric3d_z_tn_edge_probabilities(topology, p)
        direct = direct_factor_distribution(length(plan.detector_coordinates),
            topology.detector_endpoints, topology.logical_masks, probabilities)
        syndromes = BitVector[falses(length(plan.detector_coordinates))]
        for edge in plan.edges[1:4]
            syndrome = falses(length(plan.detector_coordinates))
            syndrome[collect(edge.endpoints)] .= true
            push!(syndromes, syndrome)
        end
        for syndrome in syndromes
            expected = direct_factor_weights(direct,
                length(plan.detector_coordinates), syndrome)
            actual = row_plan_weights(plan, probabilities, syndrome)
            @test actual ≈ expected atol=5e-15 rtol=5e-13
        end
    end
end

@testset "row-plan compiler rejects nonlocal repaired-transfer edges" begin
    context = build_context(2)
    topology = context.transfer_topology
    endpoints = copy(topology.detector_endpoints)
    endpoints[1] = (1, 2) # Same detector plane; impossible for this circuit.
    @test_throws ArgumentError compile_repaired_seed_transfer_plan(
        endpoints, topology.logical_masks, 2)
    endpoints[1] = reverse(topology.detector_endpoints[1])
    @test_throws ArgumentError compile_repaired_seed_transfer_plan(
        endpoints, topology.logical_masks, 2)
    endpoints = copy(topology.detector_endpoints)
    masks = copy(topology.logical_masks)
    endpoints[2] = endpoints[1]
    masks[2] = masks[1]
    @test_throws ArgumentError compile_repaired_seed_transfer_plan(
        endpoints, masks, 2)
end
