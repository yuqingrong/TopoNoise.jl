using Test, TopoNoise, LinearAlgebra

include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "dense_frontier.jl"))
include(joinpath(@__DIR__, "row_strip_mps.jl"))

using .RepairedSeedTN
using .TNDecoder3DDenseFrontier
using .TNMLRowStripMPS

BLAS.set_num_threads(1)

Threads.nthreads() > 1 || error(
    "this concurrency test must be run with julia --threads=4")

@testset "L2 Fourier sectors are deterministic across Julia threads" begin
    context = build_context(2)
    topology = context.transfer_topology
    syndrome = falses(size(topology.raw_detector_matrix, 1))
    for edge_id in (2, 7)
        for detector in topology.detector_endpoints[edge_id]
            syndrome[detector] = !syndrome[detector]
        end
    end
    exact = dense_exact_probabilities(topology, syndrome, 0.03; max_width=22)

    serial = decode_row_strip_mps(topology, 2, 0.03, syndrome;
        maxdim=4096, cutoff=0.0, algorithm=:row_mpo_fourier_sectors,
        sector_parallel=false)
    threaded_1 = decode_row_strip_mps(topology, 2, 0.03, syndrome;
        maxdim=4096, cutoff=0.0, algorithm=:row_mpo_fourier_sectors,
        sector_parallel=true)
    threaded_2 = decode_row_strip_mps(topology, 2, 0.03, syndrome;
        maxdim=4096, cutoff=0.0, algorithm=:row_mpo_fourier_sectors,
        sector_parallel=true)

    @test threaded_1.probabilities ≈ exact.probabilities atol=5e-12 rtol=5e-11
    @test isapprox(threaded_1.log_syndrome,
        exact.log_marginal_probability; atol=5e-12, rtol=5e-11)
    @test threaded_1.probabilities ≈ serial.probabilities atol=2e-14 rtol=2e-13
    @test threaded_1.log_syndrome ≈ serial.log_syndrome atol=2e-14 rtol=2e-13
    @test threaded_2.probabilities == threaded_1.probabilities
    @test threaded_2.log_syndrome == threaded_1.log_syndrome
    @test threaded_1.diagnostics.sector_parallel
    @test threaded_1.diagnostics.julia_threads == Threads.nthreads()
    @test !threaded_1.diagnostics.untruncated
    @test threaded_1.diagnostics.boundary_rank_guard_satisfied
    @test BLAS.get_num_threads() == 1
end
