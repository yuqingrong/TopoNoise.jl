using LinearAlgebra, Printf, TopoNoise, ITensors, ITensorMPS

include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "row_strip_mps.jl"))
include(joinpath(@__DIR__, "row_strip_sparse_microkernel.jl"))

using .RepairedSeedTN
using .TNMLRowStripMPS
using .TNMLRowStripSparseMicrokernel

BLAS.set_num_threads(1)

const L = 3
const P = 0.02
const SECTOR = 2
const STRIP = 9
const SEAM = 0
const MAXDIM = 240
const CUTOFF = 1e-18

function normalized_pre_row_boundary(decoder, syndrome)
    boundary = TNMLRowStripMPS._initial_sector_boundary(decoder)
    for strip_id in 1:STRIP-1
        rows = compile_fourier_fixed_seam_row_mpos(
            decoder, strip_id, syndrome, SECTOR)
        boundary = TNMLRowStripMPS._apply_fourier_fixed_seam_rows(
            decoder, boundary, rows)
        magnitude = ITensorMPS.norm(boundary.state)
        isfinite(magnitude) && magnitude > 0 || error(
            "invalid boundary norm after strip $strip_id")
        boundary.state[1] /= magnitude
    end
    return boundary
end

function main()
    context = build_context(L)
    decoder = build_row_strip_mps_decoder(context.transfer_topology, L, P;
        maxdim=MAXDIM, cutoff=CUTOFF)
    syndrome = falses(length(decoder.plan.detector_coordinates))

    println("preparing the committed ITensor boundary through strip 8")
    preparation = @timed normalized_pre_row_boundary(decoder, syndrome)
    boundary = preparation.value
    rows = compile_fourier_fixed_seam_row_mpos(
        decoder, STRIP, syndrome, SECTOR)
    row = only(filter(candidate -> candidate.seam == SEAM, rows))
    padded = TNMLRowStripMPS._pad_boundary(
        decoder, boundary, row.union_detectors, STRIP)
    operator = compile_sparse_row_operator(row, padded.state)
    input = dense_mps_chain(padded.state)
    workspace = SparseZipupWorkspace(operator, input; maxdim=MAXDIM)

    # Compilation and LAPACK initialization are outside the gate measurement.
    sparse_zipup_apply!(workspace, operator, input;
        maxdim=MAXDIM, cutoff=CUTOFF)
    GC.gc()
    measured = @timed sparse_zipup_apply!(workspace, operator, input;
        maxdim=MAXDIM, cutoff=CUTOFF)
    profile = SparseZipupTimings()
    profiled = @timed sparse_zipup_apply!(workspace, operator, input;
        maxdim=MAXDIM, cutoff=CUTOFF, timings=profile)

    println("computing the committed ITensor reference outside the timed region")
    reference_timing = @timed TNMLRowStripMPS._apply_row_operator(
        decoder, padded.state, row)
    reference = reference_timing.value
    actual_vector = dense_chain_vector(measured.value)
    reference_vector = TNMLRowStripMPS._sector_boundary_vector(
        decoder, TNMLRowStripMPS._SectorBoundaryMPS(
            reference, row.union_detectors))
    relative_error = norm(actual_vector - reference_vector) /
        norm(reference_vector)
    lognorm_error = abs(log(norm(actual_vector)) -
        log(norm(reference_vector)))

    @printf("fixture L=%d p=%.4f sector=%d strip=%d seam=%d sites=%d\n",
        L, P, SECTOR, STRIP, SEAM, length(operator.sites))
    @printf("sparse transitions=%d maximum_per_input=%d\n",
        operator.transition_count, operator.maximum_transitions_per_input)
    @printf("preparation_seconds=%.6f reference_seconds=%.6f\n",
        preparation.time, reference_timing.time)
    @printf("sparse_seconds=%.6f sparse_allocated_bytes=%d\n",
        measured.time, measured.bytes)
    @printf("profiled_repeat_seconds=%.6f profiled_allocated_bytes=%d\n",
        profiled.time, profiled.bytes)
    @printf("profile accumulation=%.6f svd=%.6f compression=%.6f final=%.6f\n",
        profile.accumulation_ns / 1e9,
        profile.factorization_ns / 1e9,
        profile.compression_ns / 1e9,
        profile.final_contraction_ns / 1e9)
    @printf("workspace_bytes=%d maximum_output_bond=%d\n",
        Base.summarysize(workspace), maximum(workspace.actual_bonds))
    @printf("relative_vector_error=%.12e lognorm_error=%.12e\n",
        relative_error, lognorm_error)
    passed_speed = measured.time <= 0.25
    passed_memory = measured.bytes <= 150_000_000
    println("speed_gate=", passed_speed,
        " memory_gate=", passed_memory,
        " overall_gate=", passed_speed && passed_memory)
    return nothing
end

main()

