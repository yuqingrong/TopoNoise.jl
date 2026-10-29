using Test, TopoNoise, ITensors, ITensorMPS

include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "dense_frontier.jl"))
include(joinpath(@__DIR__, "row_strip_mps.jl"))

using .RepairedSeedTN
using .TNDecoder3DDenseFrontier
using .TNMLRowStripPlan
using .TNMLRowStripMPS

function basis_position(detectors, bits, logical)
    mask = 0
    for (position, detector) in enumerate(detectors)
        get(bits, detector, false) && (mask |= 1 << (position - 1))
    end
    return mask + (Int(logical) << length(detectors)) + 1
end

function brute_force_row_column(plan, strip, probabilities, syndrome,
        input_mask, input_logical)
    before = [slot.detector for slot in strip.physical_before]
    after = [slot.detector for slot in strip.physical_after]
    output = zeros(Float64, 4 * (1 << length(after)))
    initial = Dict(detector => !iszero(input_mask & (1 << (position - 1)))
        for (position, detector) in enumerate(before))
    edge_ids = strip.owned_edges
    for subset in 0:(1 << length(edge_ids)) - 1
        parity = copy(initial)
        logical = UInt8(input_logical)
        weight = 1.0
        for (position, edge_id) in enumerate(edge_ids)
            edge = plan.edges[edge_id]
            present = !iszero(subset & (1 << (position - 1)))
            probability = probabilities[edge_id]
            weight *= present ? probability : 1 - probability
            if present
                for detector in edge.endpoints
                    parity[detector] = !get(parity, detector, false)
                end
                logical = xor(logical, edge.logical_mask)
            end
        end
        all(get(parity, detector, false) == Bool(syndrome[detector])
            for detector in strip.retired_detectors) || continue
        output_bits = Dict(detector => get(parity, detector, false)
            for detector in after)
        output[basis_position(after, output_bits, logical)] += weight
    end
    return output
end

@testset "L2 dense row transfer and periodic seam" begin
    context = build_context(2)
    plan = compile_repaired_seed_transfer_plan(context.transfer_topology, 2)
    probabilities = toric3d_z_tn_edge_probabilities(
        context.transfer_topology, 0.07)
    strip = plan.strips[1]
    syndrome = falses(length(plan.detector_coordinates))

    row = compile_dense_row_map(plan, strip.id, probabilities, syndrome)
    @test row.input_detectors ==
        [slot.detector for slot in strip.physical_before]
    @test row.output_detectors ==
        [slot.detector for slot in strip.physical_after]
    @test size(row.matrix) ==
        (4 * (1 << length(row.output_detectors)),
         4 * (1 << length(row.input_detectors)))
    for logical in 0:3
        expected = brute_force_row_column(plan, strip, probabilities,
            syndrome, 0, logical)
        column = basis_position(row.input_detectors, Dict{Int,Bool}(), logical)
        @test row.matrix[:, column] ≈ expected atol=2e-15 rtol=2e-14
    end

    seam_id = only(edge_id for edge_id in strip.owned_edges if begin
        edge = plan.edges[edge_id]
        source = plan.detector_coordinates[edge.endpoints[1]]
        destination = plan.detector_coordinates[edge.endpoints[2]]
        source.column == plan.L - 1 && destination.column == 0 &&
            edge.offset[3] == 1 && edge.offset[1] == 1
    end)
    seam = plan.edges[seam_id]
    deterministic = zeros(length(plan.edges))
    deterministic[seam_id] = 1.0
    seam_syndrome = falses(length(plan.detector_coordinates))
    seam_syndrome[seam.endpoints[1]] = true
    seam_row = compile_dense_row_map(
        plan, strip.id, deterministic, seam_syndrome)
    destination_bits = Dict(seam.endpoints[2] => true)
    expected_position = basis_position(seam_row.output_detectors,
        destination_bits, seam.logical_mask)
    @test seam_row.matrix[expected_position, 1] == 1.0
    @test sum(seam_row.matrix[:, 1]) == 1.0
end

@testset "L2 dense row map agrees on a nonempty boundary" begin
    context = build_context(2)
    plan = compile_repaired_seed_transfer_plan(context.transfer_topology, 2)
    probabilities = toric3d_z_tn_edge_probabilities(
        context.transfer_topology, 0.03)
    strip = plan.strips[3]
    syndrome = falses(length(plan.detector_coordinates))
    syndrome[strip.detectors[1]] = true
    row = compile_dense_row_map(plan, strip.id, probabilities, syndrome)
    masks = unique([0, 1, (1 << length(row.input_detectors)) - 1,
        0x55 & ((1 << length(row.input_detectors)) - 1)])
    for mask in masks, logical in 0:3
        expected = brute_force_row_column(plan, strip, probabilities,
            syndrome, mask, logical)
        column = mask + (logical << length(row.input_detectors)) + 1
        @test row.matrix[:, column] ≈ expected atol=3e-15 rtol=3e-14
    end
end

@testset "L2 open-chain row MPO agrees with dense row transfer" begin
    context = build_context(2)
    topology = context.transfer_topology
    decoder = build_row_strip_mps_decoder(
        topology, 2, 0.03; maxdim=4096, cutoff=0.0)
    syndrome = falses(length(decoder.plan.detector_coordinates))

    for strip_id in (1, 3, 6, 8)
        strip = decoder.plan.strips[strip_id]
        row = compile_row_mpo(decoder, strip_id, syndrome)
        @test row.union_detectors == sort(unique(vcat(
            [slot.detector for slot in strip.physical_before],
            [slot.detector for slot in strip.physical_after],
            strip.detectors));
            by=detector -> begin
                coordinate = decoder.plan.detector_coordinates[detector]
                (coordinate.column, coordinate.id)
            end)
        @test row.maximum_automaton_dimension <= 128
        input_dimension = 4 * (1 << length(strip.physical_before))
        input = collect(1.0:input_dimension)
        input ./= sum(input)
        dense = compile_dense_row_map(
            decoder.plan, strip_id, decoder.probabilities, syndrome)
        expected = dense.matrix * input
        actual = row_mpo_action(decoder, row, input)
        @test actual ≈ expected atol=3e-12 rtol=3e-11
    end
end

@testset "L2 fixed-seam row MPOs sum to the dense row transfer" begin
    context = build_context(2)
    topology = context.transfer_topology
    decoder = build_row_strip_mps_decoder(
        topology, 2, 0.03; maxdim=4096, cutoff=0.0)
    syndrome = falses(length(decoder.plan.detector_coordinates))

    for strip_id in (1, 3, 6, 8)
        strip = decoder.plan.strips[strip_id]
        rows = compile_fixed_seam_row_mpos(decoder, strip_id, syndrome)
        branch_count = 1 << length(first(rows).cross_offsets)
        @test length(rows) == branch_count
        @test getfield.(rows, :seam) ==
            Union{Nothing,Int}[0:(branch_count - 1)...]
        @test all(row -> row.maximum_automaton_dimension <= 32, rows)

        input_dimension = 4 * (1 << length(strip.physical_before))
        input = collect(1.0:input_dimension)
        input ./= sum(input)
        dense = compile_dense_row_map(
            decoder.plan, strip_id, decoder.probabilities, syndrome)
        expected = dense.matrix * input
        actual = sum(row_mpo_action(decoder, row, input) for row in rows)
        @test actual ≈ expected atol=3e-12 rtol=3e-11
    end
end

@testset "untruncated L2 row-strip MPS reproduces exact ML weights" begin
    context = build_context(2)
    topology = context.transfer_topology
    for p in (0.02, 0.10)
        syndromes = BitVector[falses(size(topology.raw_detector_matrix, 1))]
        for selection in ((1,), (2, 7), (3, 5, 11))
            syndrome = falses(size(topology.raw_detector_matrix, 1))
            for edge_id in selection
                for detector in topology.detector_endpoints[edge_id]
                    syndrome[detector] = !syndrome[detector]
                end
            end
            push!(syndromes, syndrome)
        end
        for syndrome in syndromes
            expected = dense_exact_probabilities(
                topology, syndrome, p; max_width=22)
            actual = decode_row_strip_mps(topology, 2, p, syndrome;
                maxdim=4096, cutoff=0.0, algorithm=:dense)
            @test isapprox(actual.probabilities, expected.probabilities;
                atol=2e-12, rtol=2e-11)
            @test isapprox(actual.log_syndrome,
                expected.log_marginal_probability; atol=2e-12, rtol=2e-11)
            @test actual.logical_mask == argmax(expected.probabilities) - 1
            @test actual.diagnostics.untruncated
            @test actual.diagnostics.row_updates == 8
            @test actual.diagnostics.svd_events < 8 * 16
        end
    end
end


@testset "untruncated L2 open-chain row MPO reproduces exact ML weights" begin
    context = build_context(2)
    topology = context.transfer_topology
    for p in (0.02, 0.10)
        for selection in ((), (1,), (2, 7), (3, 5, 11))
            syndrome = falses(size(topology.raw_detector_matrix, 1))
            for edge_id in selection
                for detector in topology.detector_endpoints[edge_id]
                    syndrome[detector] = !syndrome[detector]
                end
            end
            expected = dense_exact_probabilities(
                topology, syndrome, p; max_width=22)
            actual = decode_row_strip_mps(topology, 2, p, syndrome;
                maxdim=4096, cutoff=0.0, algorithm=:row_mpo)
            @test isapprox(actual.probabilities, expected.probabilities;
                atol=3e-12, rtol=3e-11)
            @test isapprox(actual.log_syndrome,
                expected.log_marginal_probability; atol=3e-12, rtol=3e-11)
            @test actual.logical_mask == argmax(expected.probabilities) - 1
            @test actual.diagnostics.algorithm ==
                "open_chain_row_mpo_boundary_mps"
            @test actual.diagnostics.row_updates == 8
            @test actual.diagnostics.mpo_applications == 8
        end
    end
end

@testset "untruncated L2 fixed-seam row MPOs reproduce exact ML weights" begin
    context = build_context(2)
    topology = context.transfer_topology
    for p in (0.02, 0.10)
        for selection in ((), (1,), (2, 7), (3, 5, 11))
            syndrome = falses(size(topology.raw_detector_matrix, 1))
            for edge_id in selection
                for detector in topology.detector_endpoints[edge_id]
                    syndrome[detector] = !syndrome[detector]
                end
            end
            expected = dense_exact_probabilities(
                topology, syndrome, p; max_width=22)
            actual = decode_row_strip_mps(topology, 2, p, syndrome;
                maxdim=4096, cutoff=0.0,
                algorithm=:row_mpo_fixed_seams)
            @test isapprox(actual.probabilities, expected.probabilities;
                atol=4e-12, rtol=4e-11)
            @test isapprox(actual.log_syndrome,
                expected.log_marginal_probability;
                atol=4e-12, rtol=4e-11)
            @test actual.logical_mask == argmax(expected.probabilities) - 1
            @test actual.diagnostics.algorithm ==
                "fixed_seam_open_chain_row_mpo_boundary_mps"
            @test actual.diagnostics.maximum_automaton_dimension <= 32
            # Bulk rows have four fixed-seam branches.  Late rows have only
            # two or one because their compiled forward stencil has shrunk.
            @test actual.diagnostics.mpo_applications == 22
        end
    end
end

@testset "L2 Fourier-sector row MPO removes the logical carrier" begin
    @test TNMLRowStripMPS._fourier_character(3, 3) == 1.0
    @test TNMLRowStripMPS._fourier_character(1, 3) == -1.0
    context = build_context(2)
    decoder = build_row_strip_mps_decoder(
        context.transfer_topology, 2, 0.03; maxdim=4096, cutoff=0.0)
    syndrome = falses(length(decoder.plan.detector_coordinates))
    for sector in 0:3, strip_id in (1, 3, 6, 8)
        rows = compile_fourier_fixed_seam_row_mpos(
            decoder, strip_id, syndrome, sector)
        @test all(row -> row.maximum_automaton_dimension <= 8, rows)
        @test all(row -> row.sector == sector, rows)
        @test all(row -> dim(last(ITensorMPS.siteinds(row.operator))) == 1,
            rows)
    end
end

@testset "untruncated L2 Fourier-sector decoder reproduces exact ML weights" begin
    context = build_context(2)
    topology = context.transfer_topology
    for p in (0.02, 0.10)
        for selection in ((), (1,), (2, 7), (3, 5, 11))
            syndrome = falses(size(topology.raw_detector_matrix, 1))
            for edge_id in selection
                for detector in topology.detector_endpoints[edge_id]
                    syndrome[detector] = !syndrome[detector]
                end
            end
            expected = dense_exact_probabilities(
                topology, syndrome, p; max_width=22)
            actual = decode_row_strip_mps(topology, 2, p, syndrome;
                maxdim=4096, cutoff=0.0,
                algorithm=:row_mpo_fourier_sectors)
            @test isapprox(actual.probabilities, expected.probabilities;
                atol=5e-12, rtol=5e-11)
            @test isapprox(actual.log_syndrome,
                expected.log_marginal_probability;
                atol=5e-12, rtol=5e-11)
            @test actual.logical_mask == argmax(expected.probabilities) - 1
            @test actual.diagnostics.algorithm ==
                "fourier_sector_fixed_seam_row_mpo_boundary_mps"
            @test actual.diagnostics.maximum_automaton_dimension <= 8
            @test actual.diagnostics.sector_count == 4
            @test actual.diagnostics.sector_signs[1] == 1
            @test actual.diagnostics.sector_logabs[1] ≈
                expected.log_marginal_probability atol=5e-12 rtol=5e-11
            @test actual.diagnostics.maximum_fourier_ratio <= 1 + 1e-10
            @test actual.diagnostics.roundoff_correction_l1 <= 1e-10
            @test sum(actual.raw_probabilities) ≈ 1.0 atol=2e-14
        end
    end
end
