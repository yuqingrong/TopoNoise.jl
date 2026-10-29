using Test, LinearAlgebra, ITensors, ITensorMPS, TopoNoise

include(joinpath(@__DIR__, "update_algorithm_benchmark.jl"))
using .TNMLUpdateAlgorithmBenchmark

function fixture(; maxdim=256)
    context = TNMLUpdateAlgorithmBenchmark.Model.build_context(2)
    settings = Toric3DZTNSettings(max_bond_dimension=maxdim, cutoff=0.0)
    decoder = TNMLUpdateAlgorithmBenchmark.Model.build_decoder(
        context, 0.02; settings).transfer
    psi = TNMLUpdateAlgorithmBenchmark.initial_boundary(decoder)
    psi = TopoNoise._toric3d_z_tn_fuse_current(psi)
    groups = decoder.groups[1]
    star = TopoNoise._toric3d_z_tn_star_mpo(psi, groups, 1, false)
    return (; psi, star)
end

function relative_tensor_distance(left, right)
    ll = real(dot(left, left))
    rr = real(dot(right, right))
    lr = real(dot(left, right))
    return sqrt(max(0.0, ll + rr - 2lr)) / max(sqrt(ll), sqrt(rr))
end

@testset "finite-prefix update algorithms" begin
    f = fixture()
    support = length(f.star)
    @test support < length(f.psi)
    boundary = linkind(f.psi, support)
    @test !isnothing(boundary)

    current, current_stats = apply_update(
        f.psi, f.star, :current; maxdim=256, cutoff=0.0)
    zipup, zipup_stats = apply_update(
        f.psi, f.star, :zipup; maxdim=256, cutoff=0.0)

    @test current_stats.svd_events == support - 1
    @test zipup_stats.svd_events == 2 * (support - 1)
    @test any(==(boundary), inds(zipup[support]))
    @test any(==(boundary), inds(zipup[support + 1]))
    @test relative_tensor_distance(current, zipup) ≤ 2e-11
    @test current_stats.elapsed_seconds ≥ 0
    @test zipup_stats.elapsed_seconds ≥ 0

    prefix = MPS(ITensor[f.psi[k] for k in 1:support];
        ortho_lims=1:support)
    package_densitymatrix = apply(f.star, prefix;
        alg="densitymatrix", maxdim=256, cutoff=0.0)
    @test !any(==(boundary), inds(package_densitymatrix[support]))

    error = try
        apply_update(f.psi, f.star, :densitymatrix; maxdim=256, cutoff=0.0)
        nothing
    catch caught
        caught
    end
    @test error isa UnsupportedFinitePrefix
    @test occursin("live suffix bond", sprint(showerror, error))
    @test_throws ArgumentError apply_update(
        f.psi, f.star, :unknown; maxdim=256, cutoff=0.0)
end

@testset "exact small decoder comparison" begin
    report = compare_small_decoder(; L=2, p=0.02, maxdim=512,
        cutoff=0.0, syndrome=falses(16))
    @test report["densitymatrix"]["status"] == "unsupported_finite_prefix"
    @test report["current"]["status"] == "complete"
    @test report["zipup"]["status"] == "complete"
    @test report["comparison"]["posterior_linf"] ≤ 2e-10
    @test report["comparison"]["abs_log_syndrome_delta"] ≤ 2e-10
    @test report["comparison"]["same_logical_mask"]
    @test report["current"]["svd_events"] > 0
    @test report["zipup"]["svd_events"] == 2report["current"]["svd_events"]
end

@testset "replay result assembly" begin
    seed = (valid=true, probabilities=[0.7, 0.1, 0.1, 0.1],
        logical_mask=0, raw_argmax=0, log_syndrome=-2.0,
        elapsed_seconds=1.25, update_seconds=1.0,
        svd_events=10, callback_events=10,
        discarded_weight=1e-9, maximum_local_truncation_error=1e-9,
        actual_maxbond=8, peak_rss_bytes=100)
    transfer = merge(seed, (probabilities=[0.6, 0.2, 0.1, 0.1],
        log_syndrome=-3.0, elapsed_seconds=2.75, svd_events=20,
        update_seconds=2.0, callback_events=20,
        actual_maxbond=16, peak_rss_bytes=200))
    result = combine_component_results(seed, transfer, 1, 1e-10)
    @test result.valid
    @test result.log_joint_syndrome == -5.0
    @test result.elapsed_seconds == 4.0
    @test result.update_seconds == 3.0
    @test result.svd_events == 30
    @test result.actual_maxbond == 16
    @test result.peak_rss_bytes == 200
    @test result.probabilities ≈ [0.22, 0.46, 0.16, 0.16]
    @test result.logical_mask == 1
end
