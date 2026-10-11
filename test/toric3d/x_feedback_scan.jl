using Test
using Random
using TopoNoise

@testset "Encoded diagonal X scan is paired, reproducible, and resumable" begin
    args = (sizes=[(2, 2, 2)], error_rates=[0.0, 0.01], shots=20,
        seed=1234, batch_size=7, routing=:one_step)
    progress = IOBuffer()
    first_scan = compare_toric3d_x_feedback_threshold(MersenneTwister(1234);
        args..., progress_io=progress)
    second_scan = compare_toric3d_x_feedback_threshold(MersenneTwister(1234);
        args...)

    @test first_scan.dimensions == [(2, 2, 2)]
    @test first_scan.error_rates == [0.0, 0.01]
    @test first_scan.shots == 20
    @test first_scan.batch_size == 7
    @test first_scan.seed == 1234
    @test first_scan.routing === :one_step
    @test length(first_scan.points) == 2
    @test length(first_scan.point_seeds) == 2
    @test first_scan.point_seeds == second_scan.point_seeds
    @test length(filter(!isempty, split(String(take!(progress)), '\n'))) == 2

    for (a, b) in zip(first_scan.points, second_scan.points)
        @test a.dimensions == b.dimensions == (2, 2, 2)
        @test a.p == b.p
        @test a.shots == b.shots == 20
        @test a.point_seed == b.point_seed
        @test a.routing === b.routing === :one_step
        @test a.metrics == b.metrics
        @test !isempty(a.decoder_version)
        for metric in values(a.metrics)
            @test metric.final_only_failures - metric.both_failures ==
                metric.final_only_only
            @test metric.layer_feedback_failures - metric.both_failures ==
                metric.layer_feedback_only
            @test 0 <= metric.both_failures <=
                min(metric.final_only_failures, metric.layer_feedback_failures) <= 20
        end
        if a.p == 0
            @test all(metric.final_only_failures == 0 &&
                metric.layer_feedback_failures == 0 for metric in values(a.metrics))
        end
    end

    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[], shots=1)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01, 0.01], shots=1)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=0)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        batch_size=0)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        routing=:invalid)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        output_dir=1)
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        output_dir="")
    @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
        MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01], shots=1,
        resume=true)

    mktempdir() do directory
        resumable = (sizes=[(2, 2, 2)], error_rates=[0.01], shots=20,
            seed=4321, batch_size=7, routing=:one_step)
        @test_throws InterruptException compare_toric3d_x_feedback_threshold(
            MersenneTwister(4321); resumable..., output_dir=directory,
            _test_stop_after_batches=1)
        checkpoint = joinpath(directory, "toric3d-x-feedback-checkpoint.toml")
        @test isfile(checkpoint)
        saved = TopoNoise.TOML.parsefile(checkpoint)
        @test saved["states"][1]["completed_shots"] == 7

        resumed = compare_toric3d_x_feedback_threshold(MersenneTwister(999);
            resumable..., output_dir=directory, resume=true)
        uninterrupted = compare_toric3d_x_feedback_threshold(MersenneTwister(4321);
            resumable...)
        @test resumed.point_seeds == uninterrupted.point_seeds
        @test [point.metrics for point in resumed.points] ==
            [point.metrics for point in uninterrupted.points]

        @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
            MersenneTwister(4321); resumable..., output_dir=directory,
            resume=true, shots=21)
        @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
            MersenneTwister(4321); resumable..., output_dir=directory,
            resume=true, routing=:full_path)
        @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
            MersenneTwister(4321); resumable..., output_dir=directory,
            resume=true, seed=9999)
        changed = Toric3DDecoderSettings(max_iter=49)
        @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
            MersenneTwister(4321); resumable..., output_dir=directory,
            resume=true, settings=changed)
    end

    mktemp() do path, io
        close(io)
        @test_throws ArgumentError compare_toric3d_x_feedback_threshold(
            MersenneTwister(1); sizes=[(2, 2, 2)], error_rates=[0.01],
            shots=1, output_dir=path)
    end
end
