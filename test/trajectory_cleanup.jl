using CairoMakie: Figure
using Random
using Test
using TopoNoise

@testset "Raw trajectory analysis cleanup" begin
    retired_symbols = (
        :MarginalSpinObservables,
        :MarginalSpinScanPoint,
        :marginal_spin_observables,
        :estimate_binder_crossings,
    )
    @test all(symbol -> !isdefined(TopoNoise, symbol), retired_symbols)

    trajectory = sample_trajectory(
        MersenneTwister(701), ToricCodeTrajectoryModel(2, 2);
        error_rate=0.25)
    @test fieldnames(TrajectoryObservables) == (
        :sampled_error_density,
        :boundary_error_density,
        :mismatch_density,
        :frustration_density,
        :largest_cluster_fraction,
        :spans_horizontal,
        :spans_vertical,
    )
    @test trajectory_observables(trajectory) isa TrajectoryObservables

    scan = scan_trajectories(
        MersenneTwister(702), [2, 3], [0.0, 1.0]; shots=8, batches=2)
    @test fieldnames(TrajectoryScanPoint) == (
        :size,
        :error_rate,
        :shots,
        :sampled_error_mean,
        :sampled_error_se,
        :boundary_error_mean,
        :boundary_error_se,
        :mismatch_mean,
        :mismatch_se,
        :frustration_mean,
        :frustration_se,
        :largest_cluster_mean,
        :largest_cluster_se,
        :horizontal_spanning_mean,
        :horizontal_spanning_se,
        :vertical_spanning_mean,
        :vertical_spanning_se,
        :horizontal_span_batches,
    )
    @test fieldnames(TrajectoryScan) == (:sizes, :error_rates, :points)
    @test_throws MethodError scan_trajectories(
        MersenneTwister(703), [2, 3], [0.0, 1.0];
        shots=8, batches=2, spin_samples=1)

    crossings = estimate_crossings(MersenneTwister(704), scan; bootstrap=4)
    figure = plot_trajectory_scan(scan, crossings)
    @test figure isa Figure
    titles = [string(block.title[]) for block in figure.content
              if hasproperty(block, :title)]
    @test Set(titles) == Set((
        "Injected errors and mismatches",
        "Plaquette frustration",
        "Largest occupied cluster",
        "Horizontal spanning probability",
    ))
end
