using Random
using Test
using TopoNoise

@testset "Welford streaming accumulator" begin
    accumulator = TopoNoise._WelfordAccumulator()
    for value in (1.0, 2.0, 3.0, 4.0)
        TopoNoise._push!(accumulator, value)
    end
    @test accumulator.count == 4
    @test accumulator.mean == 2.5
    @test TopoNoise._standard_error(accumulator) ≈ sqrt(5 / 12)
end

@testset "Streaming raw trajectory scans" begin
    scan = scan_trajectories(
        MersenneTwister(201), [2, 3], [0.0, 1.0]; shots=20, batches=4)
    @test scan isa TrajectoryScan
    @test scan.sizes == [2, 3]
    @test scan.error_rates == [0.0, 1.0]
    @test length(scan.points) == 4

    for size in scan.sizes
        zero_point = only(filter(
            point -> point.size == size && point.error_rate == 0.0,
            scan.points))
        @test zero_point.shots == 20
        @test zero_point.sampled_error_mean == 0.0
        @test zero_point.boundary_error_mean == 0.0
        @test zero_point.mismatch_mean == 0.0
        @test zero_point.frustration_mean == 0.0
        @test zero_point.largest_cluster_mean ≈ 1 / size^2
        @test zero_point.horizontal_spanning_mean == 0.0
        @test zero_point.vertical_spanning_mean == 0.0
        @test zero_point.horizontal_span_batches == zeros(4)

        one_point = only(filter(
            point -> point.size == size && point.error_rate == 1.0,
            scan.points))
        @test one_point.sampled_error_mean == 1.0
        @test one_point.boundary_error_mean == 1.0
        @test one_point.mismatch_mean == 1.0
        @test one_point.frustration_mean == 0.0
        @test one_point.largest_cluster_mean == 1.0
        @test one_point.horizontal_spanning_mean == 1.0
        @test one_point.vertical_spanning_mean == 1.0
        @test one_point.horizontal_span_batches == ones(4)
    end

    repeated = scan_trajectories(
        MersenneTwister(201), [2, 3], [0.0, 1.0]; shots=20, batches=4)
    @test repeated.sizes == scan.sizes
    @test repeated.error_rates == scan.error_rates
    for (first, second) in zip(repeated.points, scan.points)
        for field in fieldnames(TrajectoryScanPoint)
            @test getfield(first, field) == getfield(second, field)
        end
    end

    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [1, 2], [0.5]; shots=10, batches=2)
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [3, 2], [0.5]; shots=10, batches=2)
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [2, 2], [0.5]; shots=10, batches=2)
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [2], [0.6, 0.5]; shots=10, batches=2)
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [2], [0.5]; shots=0, batches=2)
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [2], [0.5]; shots=10, batches=11)
end

function synthetic_scan_point(size, error_rate, horizontal_probability)
    batches = fill(Float64(horizontal_probability), 4)
    return TrajectoryScanPoint(
        size, Float64(error_rate), 100,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        Float64(horizontal_probability), 0.0,
        0.0, 0.0,
        batches)
end

function with_horizontal_batches(point, batches)
    return TrajectoryScanPoint(
        point.size, point.error_rate, point.shots,
        point.sampled_error_mean, point.sampled_error_se,
        point.boundary_error_mean, point.boundary_error_se,
        point.mismatch_mean, point.mismatch_se,
        point.frustration_mean, point.frustration_se,
        point.largest_cluster_mean, point.largest_cluster_se,
        point.horizontal_spanning_mean, point.horizontal_spanning_se,
        point.vertical_spanning_mean, point.vertical_spanning_se,
        Float64[batches...])
end

@testset "Finite-size horizontal-spanning crossings" begin
    rates = [0.4, 0.5, 0.6]
    points = TrajectoryScanPoint[
        synthetic_scan_point(4, 0.4, 0.1),
        synthetic_scan_point(4, 0.5, 0.5),
        synthetic_scan_point(4, 0.6, 0.9),
        synthetic_scan_point(8, 0.4, 0.2),
        synthetic_scan_point(8, 0.5, 0.5),
        synthetic_scan_point(8, 0.6, 0.8),
    ]
    scan = TrajectoryScan([4, 8], rates, points)
    crossing = only(estimate_crossings(
        MersenneTwister(202), scan; bootstrap=100, confidence=0.95))
    @test crossing.small_size == 4
    @test crossing.large_size == 8
    @test crossing.estimate == 0.5
    @test crossing.ci_low == 0.5
    @test crossing.ci_high == 0.5
    @test crossing.valid_bootstrap_fraction == 1.0
    @test crossing.status == :ok

    unbracketed = TrajectoryScan(
        [4, 8], [0.4, 0.6], [
            synthetic_scan_point(4, 0.4, 0.1),
            synthetic_scan_point(4, 0.6, 0.2),
            synthetic_scan_point(8, 0.4, 0.3),
            synthetic_scan_point(8, 0.6, 0.4),
        ])
    @test only(estimate_crossings(
        MersenneTwister(203), unbracketed; bootstrap=20)).status == :unbracketed

    identical = TrajectoryScan(
        [4, 8], rates, [
            synthetic_scan_point(4, 0.4, 0.1),
            synthetic_scan_point(4, 0.5, 0.5),
            synthetic_scan_point(4, 0.6, 0.9),
            synthetic_scan_point(8, 0.4, 0.1),
            synthetic_scan_point(8, 0.5, 0.5),
            synthetic_scan_point(8, 0.6, 0.9),
        ])
    @test only(estimate_crossings(
        MersenneTwister(204), identical; bootstrap=20)).status == :unstable

    one_batch_points = [with_horizontal_batches(
        point, [point.horizontal_spanning_mean]) for point in points]
    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(205),
        TrajectoryScan([4, 8], rates, one_batch_points); bootstrap=20)

    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(1), scan; bootstrap=0)
    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(1), scan; confidence=1.0)
end
