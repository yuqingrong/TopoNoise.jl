using Test
using Random
using TopoNoise

@testset "Eight joint patterns determine all three protocol comparisons" begin
    counts = [2, 1, 3, 0, 2, 1, 0, 1]
    metric = TopoNoise._toric3d_triple_metric(10, counts)
    @test metric isa Toric3DTripleMetricEstimate
    @test metric.shots == 10
    @test collect(metric.joint_counts) == counts

    @test metric.none.failures == 3
    @test metric.slice_local.failures == 4
    @test metric.full_prefix.failures == 4
    @test metric.none.rate == 0.3
    @test metric.slice_local.rate == metric.full_prefix.rate == 0.4
    @test metric.none.ci == TopoNoise._toric3d_wilson(3, 10)

    pair = metric.none_vs_slice_local
    @test (pair.neither, pair.first_only, pair.second_only, pair.both) == (4, 2, 3, 1)
    @test pair.first === :none
    @test pair.second === :slice_local
    @test pair.difference == 0.1
    @test pair.standard_error ≈ sqrt(0.49 / 9)
    @test pair.ci[1] < pair.difference < pair.ci[2]

    pair = metric.none_vs_full_prefix
    @test (pair.neither, pair.first_only, pair.second_only, pair.both) == (5, 1, 2, 2)
    @test pair.difference == 0.1
    pair = metric.slice_local_vs_full_prefix
    @test (pair.neither, pair.first_only, pair.second_only, pair.both) == (3, 3, 3, 1)
    @test pair.difference == 0.0

    @test_throws ArgumentError TopoNoise._toric3d_triple_metric(
        10, [2, 1, 3, 0, 2, 1, 0, -1])
    @test_throws ArgumentError TopoNoise._toric3d_triple_metric(
        11, counts)
    @test_throws ArgumentError TopoNoise._toric3d_triple_metric(
        10, counts[1:7])
end

@testset "Diagonal feedback scan is paired, seeded, and channel complete" begin
    args = (sizes=[(2, 2, 2)], error_rates=[0.0, 0.01],
        channels=[:x_only, :z_only], shots=20, seed=1234,
        batch_size=7, x_routing=:one_step)
    progress = IOBuffer()
    first_scan = compare_toric3d_diagonal_feedback(Xoshiro(1234);
        args..., progress_io=progress)
    second_scan = compare_toric3d_diagonal_feedback(Xoshiro(1234); args...)

    @test first_scan isa Toric3DDiagonalFeedbackScan
    @test first_scan.dimensions == [(2, 2, 2)]
    @test first_scan.error_rates == [0.0, 0.01]
    @test first_scan.channels == [:x_only, :z_only]
    @test first_scan.shots == 20
    @test first_scan.seed == 1234
    @test first_scan.batch_size == 7
    @test first_scan.x_routing === :one_step
    @test first_scan.decoder_builds == (x=2, z=1)
    @test length(first_scan.points) == 4
    @test length(first_scan.point_seeds) == 4
    @test first_scan.point_seeds == second_scan.point_seeds
    @test length(filter(!isempty, split(String(take!(progress)), '\n'))) == 4

    for (first_point, second_point) in zip(first_scan.points, second_scan.points)
        @test first_point.dimensions == second_point.dimensions == (2, 2, 2)
        @test first_point.p == second_point.p
        @test first_point.channel == second_point.channel
        @test first_point.point_seed == second_point.point_seed
        @test first_point.invalid_recoveries == 0
        @test keys(first_point.metrics) == (:logical_1, :logical_2, :state_failure)
        @test first_point.metrics == second_point.metrics
        @test !isempty(first_point.decoder_versions.pymatching)
        @test all(value >= 0 for value in values(first_point.phase_seconds))
        for metric in values(first_point.metrics)
            @test sum(metric.joint_counts) == 20
            @test metric.none.failures == sum(metric.joint_counts[index]
                for index in 1:8 if isodd(index - 1))
            @test metric.slice_local.failures == sum(metric.joint_counts[index]
                for index in 1:8 if ((index - 1) & 2) != 0)
            @test metric.full_prefix.failures == sum(metric.joint_counts[index]
                for index in 1:8 if ((index - 1) & 4) != 0)
        end
        if first_point.p == 0
            @test all(all(protocol.failures == 0 for protocol in
                (metric.none, metric.slice_local, metric.full_prefix))
                for metric in values(first_point.metrics))
        end
    end

    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[], error_rates=[0.01], shots=1)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[], shots=1)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[0.01], channels=[], shots=1)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[0.01], channels=[:bad], shots=1)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=0)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1, batch_size=0)
    @test_throws ArgumentError compare_toric3d_diagonal_feedback(
        Xoshiro(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        x_routing=:bad)
end

@testset "Interrupted diagonal feedback scans resume exact point substreams" begin
    mktempdir() do directory
        args = (sizes=[(2, 2, 2)], error_rates=[0.01],
            channels=[:x_only, :z_only], shots=20, seed=4321,
            batch_size=7, x_routing=:one_step)
        @test_throws InterruptException compare_toric3d_diagonal_feedback(
            Xoshiro(4321); args..., output_dir=directory,
            _test_stop_after_batches=1)
        checkpoint = joinpath(directory,
            "toric3d-diagonal-feedback-checkpoint.toml")
        @test isfile(checkpoint)
        saved = TopoNoise.TOML.parsefile(checkpoint)
        @test saved["states"][1]["completed_shots"] == 7

        resumed = compare_toric3d_diagonal_feedback(Xoshiro(999);
            args..., output_dir=directory, resume=true)
        uninterrupted = compare_toric3d_diagonal_feedback(Xoshiro(4321); args...)
        @test resumed.point_seeds == uninterrupted.point_seeds
        @test [point.metrics for point in resumed.points] ==
              [point.metrics for point in uninterrupted.points]

        mismatches = (
            (sizes=[(1, 2, 2)],),
            (error_rates=[0.02],),
            (channels=[:x_only],),
            (shots=21,),
            (batch_size=5,),
            (x_routing=:full_path,),
            (seed=9999,),
            (settings=Toric3DDecoderSettings(max_iter=49),),
        )
        for changed in mismatches
            draws = Ref(0)
            changed_args = merge(args, changed)
            @test_throws ArgumentError compare_toric3d_diagonal_feedback(
                Xoshiro(5); changed_args..., output_dir=directory, resume=true,
                _test_fault_draw_hook=() -> (draws[] += 1))
            @test draws[] == 0
        end
    end
end

@testset "Invalid recovery context is saved before scan abort" begin
    mktempdir() do directory
        hook = context -> falses(context.recovery_length)
        @test_throws ErrorException compare_toric3d_diagonal_feedback(
            Xoshiro(77); sizes=[(2, 2, 2)], error_rates=[0.01],
            channels=[:x_only], shots=2, seed=77, batch_size=2,
            output_dir=directory, _test_recovery_hook=hook)
        failure_path = joinpath(directory,
            "toric3d-diagonal-feedback-failure.toml")
        checkpoint_path = joinpath(directory,
            "toric3d-diagonal-feedback-checkpoint.toml")
        @test isfile(failure_path)
        failure = TopoNoise.TOML.parsefile(failure_path)
        @test failure["dimensions"] == [2, 2, 2]
        @test failure["channel"] == "x_only"
        @test failure["p"] == 0.01
        @test failure["protocol"] == "full_prefix"
        @test failure["checkpoint"] == "final"
        @test failure["batch_start"] == 1
        @test !isempty(failure["point_seed"])
        @test contains(failure["error"], "invalid X recovery")
        checkpoint = TopoNoise.TOML.parsefile(checkpoint_path)
        @test checkpoint["invalid_recoveries"] == 1
    end

    mktempdir() do directory
        compare_toric3d_diagonal_feedback(Xoshiro(79);
            sizes=[(1, 2, 2)], error_rates=[0.0], channels=[:z_only],
            shots=2, seed=79, batch_size=1, output_dir=directory)
        checkpoint = TopoNoise.TOML.parsefile(joinpath(directory,
            "toric3d-diagonal-feedback-checkpoint.toml"))
        @test checkpoint["invalid_recoveries"] == 0
        @test all(state["invalid_recoveries"] == 0 for state in checkpoint["states"])
    end
end

@testset "User interruption is not recorded as an invalid recovery" begin
    mktempdir() do directory
        hook = _ -> throw(InterruptException())
        @test_throws InterruptException compare_toric3d_diagonal_feedback(
            Xoshiro(81); sizes=[(2, 2, 2)], error_rates=[0.01],
            channels=[:x_only], shots=2, seed=81, batch_size=2,
            output_dir=directory, _test_recovery_hook=hook)
        checkpoint_path = joinpath(directory,
            "toric3d-diagonal-feedback-checkpoint.toml")
        failure_path = joinpath(directory,
            "toric3d-diagonal-feedback-failure.toml")
        @test isfile(checkpoint_path)
        checkpoint = TopoNoise.TOML.parsefile(checkpoint_path)
        @test checkpoint["invalid_recoveries"] == 0
        @test all(state["invalid_recoveries"] == 0
            for state in checkpoint["states"])
        @test !isfile(failure_path)
    end
end
