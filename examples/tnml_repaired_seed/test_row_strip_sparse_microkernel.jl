using Test, LinearAlgebra, TopoNoise, ITensors, ITensorMPS

include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "row_strip_mps.jl"))
include(joinpath(@__DIR__, "row_strip_sparse_microkernel.jl"))

using .RepairedSeedTN
using .TNMLRowStripMPS
using .TNMLRowStripSparseMicrokernel

function branch_fixture(decoder, syndrome, sector, seam)
    rows = compile_fourier_fixed_seam_row_mpos(
        decoder, 1, syndrome, sector)
    row = only(filter(candidate -> candidate.seam == seam, rows))
    boundary = TNMLRowStripMPS._initial_sector_boundary(decoder)
    padded = TNMLRowStripMPS._pad_boundary(
        decoder, boundary, row.union_detectors, row.strip_id)
    return row, padded
end

function exact_branch_vector(decoder, row, padded; maxdim=4096)
    reference = ITensorMPS.apply(row.operator, padded.state;
        alg="zipup", maxdim=maxdim, cutoff=0.0,
        truncate_kwargs=(; maxdim, cutoff=0.0))
    return TNMLRowStripMPS._sector_boundary_vector(
        decoder, TNMLRowStripMPS._SectorBoundaryMPS(
            reference, row.union_detectors))
end

function exact_schmidt_rank_bound(physical_dimensions)
    total = prod(physical_dimensions)
    left = 1
    maximum_rank = 1
    for position in 1:length(physical_dimensions)-1
        left *= physical_dimensions[position]
        maximum_rank = max(maximum_rank, min(left, total ÷ left))
    end
    return maximum_rank
end

@testset "retained-rank cutoff policy" begin
    spectrum = [4.0, 3.0, 2.0, 1.0]
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        spectrum, 2, 0.0) == 2
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        spectrum, 10, 0.0) == 4
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        spectrum, 4, 0.04) == 3
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        spectrum, 4, 0.20) == 2
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        zeros(4), 4, 0.99) == 1
    # The discarded-weight comparison includes equality at the cutoff.
    @test TNMLRowStripSparseMicrokernel._retained_rank(
        [3.0, 1.0], 2, 0.10) == 1
end

@testset "sparse row zip-up reproduces exact L2 Fourier seam branches" begin
    context = build_context(2)
    decoder = build_row_strip_mps_decoder(
        context.transfer_topology, 2, 0.02; maxdim=4096, cutoff=0.0)
    zero_syndrome = falses(length(decoder.plan.detector_coordinates))
    nonzero_syndrome = copy(zero_syndrome)
    edge = decoder.plan.edges[first(decoder.plan.strips[1].owned_edges)]
    for detector in edge.endpoints
        nonzero_syndrome[detector] = !nonzero_syndrome[detector]
    end
    @test any(nonzero_syndrome)

    # The second case simultaneously exercises a nonzero observed syndrome,
    # a different Fourier character, and a nonzero periodic seam assignment.
    for (syndrome, sector, seam) in
            ((zero_syndrome, 2, 0), (nonzero_syndrome, 1, 3))
        row, padded = branch_fixture(decoder, syndrome, sector, seam)
        operator = compile_sparse_row_operator(row, padded.state)
        input = dense_mps_chain(padded.state)
        workspace = SparseZipupWorkspace(operator, input; maxdim=4096)

        # These are algebraically exact tests: cutoff is zero, and the
        # physical-state Schmidt bound and recursively computed zip-up bounds
        # are strictly below the requested maximum. Hence neither reference
        # nor custom maxdim can truncate any local SVD.
        physical_dimensions = vcat(fill(2, length(row.union_detectors)), 1)
        @test exact_schmidt_rank_bound(physical_dimensions) < 4096
        @test all(bound -> bound < workspace.maxdim,
            workspace.bond_upper_bounds)
        actual = sparse_zipup_apply!(workspace, operator, input;
            maxdim=4096, cutoff=0.0)
        first_vector = copy(dense_chain_vector(actual))
        @test maximum(workspace.actual_bonds) < workspace.maxdim
        expected_vector = exact_branch_vector(decoder, row, padded)

        @test operator.maximum_transitions_per_input <= 8
        @test norm(expected_vector) > 0
        @test first_vector ≈ expected_vector atol=2e-12 rtol=2e-11

        # A workspace is reusable without retaining data from the previous
        # application, including for signed Fourier-sector amplitudes.
        repeated = sparse_zipup_apply!(workspace, operator, input;
            maxdim=4096, cutoff=0.0)
        @test isapprox(dense_chain_vector(repeated), first_vector;
            atol=2e-14, rtol=2e-13)
    end
end

@testset "finite maxdim and cutoff remain bounded and deterministic" begin
    context = build_context(2)
    decoder = build_row_strip_mps_decoder(
        context.transfer_topology, 2, 0.07; maxdim=2, cutoff=1e-6)
    syndrome = falses(length(decoder.plan.detector_coordinates))
    edge = decoder.plan.edges[first(decoder.plan.strips[1].owned_edges)]
    for detector in edge.endpoints
        syndrome[detector] = !syndrome[detector]
    end
    row, padded = branch_fixture(decoder, syndrome, 3, 1)
    operator = compile_sparse_row_operator(row, padded.state)
    input = dense_mps_chain(padded.state)
    workspace = SparseZipupWorkspace(operator, input; maxdim=2)

    first_result = sparse_zipup_apply!(workspace, operator, input;
        maxdim=2, cutoff=1e-6)
    first_vector = copy(dense_chain_vector(first_result))
    @test all(bond -> bond <= 2, workspace.actual_bonds)
    @test any(==(2), workspace.actual_bonds)
    @test all(isfinite, first_vector)
    @test norm(first_vector) > 0

    second_result = sparse_zipup_apply!(workspace, operator, input;
        maxdim=2, cutoff=1e-6)
    @test isapprox(dense_chain_vector(second_result), first_vector;
        atol=2e-14, rtol=2e-13)
end

