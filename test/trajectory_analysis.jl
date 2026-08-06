using Random

@testset "Welford streaming accumulator" begin
    accumulator = TopoNoise._WelfordAccumulator()
    for value in (1.0, 2.0, 3.0, 4.0)
        TopoNoise._push!(accumulator, value)
    end
    @test accumulator.count == 4
    @test accumulator.mean == 2.5
    @test TopoNoise._standard_error(accumulator) ≈ sqrt(5 / 3 / 4)
end

@testset "Streaming trajectory scans" begin
    scan = scan_trajectories(
        MersenneTwister(201), [2, 3], [0.0, 1.0]; shots=20, batches=4)
    @test scan isa TrajectoryScan
    @test scan.sizes == [2, 3]
    @test scan.error_rates == [0.0, 1.0]
    @test length(scan.points) == 4
    @test length(scan.marginal_spin_points) == 4

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

        zero_spin = only(filter(
            point -> point.size == size && point.error_rate == 0.0,
            scan.marginal_spin_points))
        @test zero_spin.absolute_magnetization_mean == 1.0
        @test zero_spin.second_moment_mean == 1.0
        @test zero_spin.fourth_moment_mean == 1.0
        @test zero_spin.binder_cumulant ≈ 2 / 3
        @test zero_spin.binder_cumulant_se == 0.0
        @test zero_spin.batch_counts == fill(5, 4)
        @test zero_spin.second_moment_batches == ones(4)
        @test zero_spin.fourth_moment_batches == ones(4)

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

        sites = size^2
        one_spin = only(filter(
            point -> point.size == size && point.error_rate == 1.0,
            scan.marginal_spin_points))
        @test one_spin.second_moment_mean == 1 / sites
        @test one_spin.fourth_moment_mean ==
              (3 * sites^2 - 2 * sites) / sites^4
        @test one_spin.binder_cumulant ≈ 2 / (3 * sites)
        @test one_spin.binder_cumulant_se ≈ 0.0 atol=1e-14
        @test one_spin.second_moment_batches ≈ fill(1 / sites, 4)
        @test one_spin.fourth_moment_batches ≈
              fill((3 * sites^2 - 2 * sites) / sites^4, 4)
    end

    repeated = scan_trajectories(
        MersenneTwister(201), [2, 3], [0.0, 1.0]; shots=20, batches=4)
    @test [point.horizontal_span_batches for point in repeated.points] ==
          [point.horizontal_span_batches for point in scan.points]

    more_spin_samples = scan_trajectories(
        MersenneTwister(201), [2, 3], [0.0, 1.0];
        shots=20, batches=4, spin_samples=7)
    for (first, second) in zip(more_spin_samples.points, scan.points)
        for field in fieldnames(TrajectoryScanPoint)
            @test getfield(first, field) == getfield(second, field)
        end
    end

    one_batch = scan_trajectories(
        MersenneTwister(209), [2], [0.5]; shots=5, batches=1)
    one_batch_spin = only(one_batch.marginal_spin_points)
    @test one_batch_spin.binder_cumulant ≈ TopoNoise._binder_cumulant(
        one_batch_spin.second_moment_mean, one_batch_spin.fourth_moment_mean)
    @test ismissing(one_batch_spin.binder_cumulant_se)

    unequal_counts = [1, 2, 3]
    unequal_second = [0.2, 0.4, 0.8]
    unequal_fourth = [0.05, 0.2, 1.2]
    jackknife_replicates = [0.34895833333333337,
                            0.28007889546351083, 0.55]
    jackknife_center = sum(jackknife_replicates) / 3
    expected_jackknife = sqrt(2 / 3 * sum(
        (replicate - jackknife_center)^2
        for replicate in jackknife_replicates))
    @test TopoNoise._binder_jackknife_se(
        unequal_counts, unequal_second, unequal_fourth) ≈
          expected_jackknife atol=1e-14

    aggregate_binder = TopoNoise._binder_cumulant(
        (1.0 + 1 / 4) / 2, (1.0 + 5 / 32) / 2)
    per_trajectory_binder = (
        TopoNoise._binder_cumulant(1.0, 1.0) +
        TopoNoise._binder_cumulant(1 / 4, 5 / 32)) / 2
    @test aggregate_binder ≈ 38 / 75
    @test per_trajectory_binder ≈ 5 / 12

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
    @test_throws ArgumentError scan_trajectories(
        MersenneTwister(1), [2], [0.5];
        shots=10, batches=2, spin_samples=0)
end

function synthetic_scan_point(size, error_rate, horizontal_probability)
    batches = fill(Float64(horizontal_probability), 4)
    return TrajectoryScanPoint(
        size, Float64(error_rate), 100,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0,
        Float64(horizontal_probability), 0.0,
        0.0, 0.0,
        batches)
end

function synthetic_spin_point(size, error_rate, binder; batches=4)
    second = 0.5
    fourth = 3 * second^2 * (1 - binder)
    return MarginalSpinScanPoint(
        size, Float64(error_rate), 100,
        0.5, 0.0,
        second, 0.0,
        fourth, 0.0,
        Float64(binder), 0.0,
        fill(div(100, batches), batches),
        fill(second, batches), fill(fourth, batches))
end

function variable_spin_point(
        size, error_rate, binder_batches; batch_counts=fill(25, length(binder_batches)))
    second_batches = fill(0.5, length(binder_batches))
    fourth_batches = 0.75 .* (1 .- binder_batches)
    shots = sum(batch_counts)
    second = sum(batch_counts .* second_batches) / shots
    fourth = sum(batch_counts .* fourth_batches) / shots
    binder = TopoNoise._binder_cumulant(second, fourth)
    return MarginalSpinScanPoint(
        size, Float64(error_rate), shots,
        0.5, 0.0,
        second, 0.0,
        fourth, 0.0,
        binder, TopoNoise._binder_jackknife_se(
            batch_counts, second_batches, fourth_batches),
        Int[batch_counts...], second_batches, fourth_batches)
end

function binder_scan(rates, small_curve, large_curve; batches=4)
    raw_points = TrajectoryScanPoint[]
    spin_points = MarginalSpinScanPoint[]
    for (size, curve) in ((4, small_curve), (8, large_curve))
        for (rate, binder) in zip(rates, curve)
            push!(raw_points, synthetic_scan_point(size, rate, 0.5))
            push!(spin_points, synthetic_spin_point(
                size, rate, binder; batches=batches))
        end
    end
    return TrajectoryScan([4, 8], Float64[rates...], raw_points, spin_points)
end

@testset "Marginalized Binder crossings" begin
    rates = [0.4, 0.5, 0.6]
    scan = binder_scan(rates, [0.60, 0.45, 0.20], [0.55, 0.45, 0.25])
    crossing = only(estimate_binder_crossings(
        MersenneTwister(301), scan; bootstrap=100, confidence=0.95))
    @test crossing.small_size == 4
    @test crossing.large_size == 8
    @test crossing.estimate == 0.5
    @test crossing.ci_low == 0.5
    @test crossing.ci_high == 0.5
    @test crossing.valid_bootstrap_fraction == 1.0
    @test crossing.status == :ok

    paired_point = MarginalSpinScanPoint(
        4, 0.5, 6,
        0.5, 0.0,
        3.4 / 6, 0.0,
        4.05 / 6, 0.0,
        TopoNoise._binder_cumulant(3.4 / 6, 4.05 / 6), 0.0,
        [1, 2, 3], [0.2, 0.4, 0.8], [0.05, 0.2, 1.2])
    expected_rng = MersenneTwister(310)
    sampled_indices = rand(expected_rng, 1:3, 3)
    sampled_count = sum(paired_point.batch_counts[sampled_indices])
    expected_second = sum(
        paired_point.batch_counts[index] *
        paired_point.second_moment_batches[index]
        for index in sampled_indices) / sampled_count
    expected_fourth = sum(
        paired_point.batch_counts[index] *
        paired_point.fourth_moment_batches[index]
        for index in sampled_indices) / sampled_count
    @test TopoNoise._bootstrap_binder(
        MersenneTwister(310), paired_point) ≈
          TopoNoise._binder_cumulant(expected_second, expected_fourth)

    noisy_points = MarginalSpinScanPoint[]
    for (size, curve) in ((4, [0.60, 0.45, 0.20]),
                          (8, [0.55, 0.45, 0.25]))
        for (rate, binder) in zip(rates, curve)
            push!(noisy_points, variable_spin_point(
                size, rate, [binder - 0.30, binder + 0.30]))
        end
    end
    noisy_scan = TrajectoryScan(
        [4, 8], Float64[rates...], scan.points, noisy_points)
    low_valid = only(estimate_binder_crossings(
        MersenneTwister(311), noisy_scan; bootstrap=400))
    @test low_valid.status == :unstable
    @test low_valid.valid_bootstrap_fraction < 0.8
    @test ismissing(low_valid.ci_low) && ismissing(low_valid.ci_high)

    unbracketed = only(estimate_binder_crossings(
        MersenneTwister(302),
        binder_scan(rates, [0.60, 0.50, 0.40], [0.50, 0.40, 0.30]);
        bootstrap=20))
    @test ismissing(unbracketed.estimate)
    @test unbracketed.status == :unbracketed

    identical = only(estimate_binder_crossings(
        MersenneTwister(303),
        binder_scan(rates, [0.60, 0.45, 0.20], [0.60, 0.45, 0.20]);
        bootstrap=20))
    @test ismissing(identical.estimate)
    @test identical.status == :unstable

    plateau_rates = [0.3, 0.4, 0.5, 0.6]
    plateau = only(estimate_binder_crossings(
        MersenneTwister(304), binder_scan(
            plateau_rates,
            [0.65, 0.50, 0.40, 0.20],
            [0.60, 0.50, 0.40, 0.25]); bootstrap=20))
    @test ismissing(plateau.estimate)
    @test plateau.status == :unstable

    multiple_rates = [0.3, 0.4, 0.5, 0.6]
    multiple = only(estimate_binder_crossings(
        MersenneTwister(305), binder_scan(
            multiple_rates,
            [0.65, 0.45, 0.35, 0.15],
            [0.60, 0.50, 0.30, 0.20]); bootstrap=20))
    @test ismissing(multiple.estimate)
    @test multiple.status == :unstable

    @test_throws ArgumentError estimate_binder_crossings(
        MersenneTwister(306),
        binder_scan(rates, [0.60, 0.45, 0.20], [0.55, 0.45, 0.25];
                    batches=1);
        bootstrap=20)
    @test_throws ArgumentError estimate_binder_crossings(
        MersenneTwister(307), scan; bootstrap=0)
    @test_throws ArgumentError estimate_binder_crossings(
        MersenneTwister(308), scan; confidence=1.0)

    raw_only = TrajectoryScan(scan.sizes, scan.error_rates, scan.points)
    @test isempty(raw_only.marginal_spin_points)
    @test_throws ArgumentError estimate_binder_crossings(
        MersenneTwister(309), raw_only; bootstrap=20)
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

@testset "Finite-size spanning crossings" begin
    rates = [0.4, 0.5, 0.6]
    points = TrajectoryScanPoint[]
    append!(points, [
        synthetic_scan_point(4, 0.4, 0.1),
        synthetic_scan_point(4, 0.5, 0.5),
        synthetic_scan_point(4, 0.6, 0.9),
        synthetic_scan_point(8, 0.4, 0.2),
        synthetic_scan_point(8, 0.5, 0.5),
        synthetic_scan_point(8, 0.6, 0.8),
    ])
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
    missing_crossing = only(estimate_crossings(
        MersenneTwister(203), unbracketed; bootstrap=20))
    @test ismissing(missing_crossing.estimate)
    @test ismissing(missing_crossing.ci_low)
    @test ismissing(missing_crossing.ci_high)
    @test missing_crossing.valid_bootstrap_fraction == 0.0
    @test missing_crossing.status == :unbracketed

    endpoint_only = TrajectoryScan(
        [4, 8], [0.0, 0.5, 1.0], [
            synthetic_scan_point(4, 0.0, 0.0),
            synthetic_scan_point(4, 0.5, 0.3),
            synthetic_scan_point(4, 1.0, 0.6),
            synthetic_scan_point(8, 0.0, 0.0),
            synthetic_scan_point(8, 0.5, 0.4),
            synthetic_scan_point(8, 1.0, 0.7),
        ])
    endpoint_crossing = only(estimate_crossings(
        MersenneTwister(204), endpoint_only; bootstrap=20))
    @test ismissing(endpoint_crossing.estimate)
    @test endpoint_crossing.status == :unbracketed

    identical = TrajectoryScan(
        [4, 8], rates, [
            synthetic_scan_point(4, 0.4, 0.1),
            synthetic_scan_point(4, 0.5, 0.5),
            synthetic_scan_point(4, 0.6, 0.9),
            synthetic_scan_point(8, 0.4, 0.1),
            synthetic_scan_point(8, 0.5, 0.5),
            synthetic_scan_point(8, 0.6, 0.9),
        ])
    identical_crossing = only(estimate_crossings(
        MersenneTwister(205), identical; bootstrap=20))
    @test ismissing(identical_crossing.estimate)
    @test identical_crossing.status == :unstable

    plateau_rates = [0.3, 0.4, 0.5, 0.6]
    plateau = TrajectoryScan(
        [4, 8], plateau_rates, [
            synthetic_scan_point(4, 0.3, 0.1),
            synthetic_scan_point(4, 0.4, 0.5),
            synthetic_scan_point(4, 0.5, 0.5),
            synthetic_scan_point(4, 0.6, 0.9),
            synthetic_scan_point(8, 0.3, 0.2),
            synthetic_scan_point(8, 0.4, 0.5),
            synthetic_scan_point(8, 0.5, 0.5),
            synthetic_scan_point(8, 0.6, 0.8),
        ])
    plateau_crossing = only(estimate_crossings(
        MersenneTwister(206), plateau; bootstrap=20))
    @test ismissing(plateau_crossing.estimate)
    @test plateau_crossing.status == :unstable

    one_batch_points = [with_horizontal_batches(
        point, [point.horizontal_spanning_mean]) for point in points]
    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(207),
        TrajectoryScan([4, 8], rates, one_batch_points);
        bootstrap=20)

    invalid_bootstrap_points = [
        with_horizontal_batches(point, fill(point.size == 4 ? 0.0 : 1.0, 4))
        for point in points
    ]
    unstable_bootstrap = only(estimate_crossings(
        MersenneTwister(208),
        TrajectoryScan([4, 8], rates, invalid_bootstrap_points);
        bootstrap=20))
    @test unstable_bootstrap.estimate == 0.5
    @test ismissing(unstable_bootstrap.ci_low)
    @test ismissing(unstable_bootstrap.ci_high)
    @test unstable_bootstrap.valid_bootstrap_fraction == 0.0
    @test unstable_bootstrap.status == :unstable

    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(1), scan; bootstrap=0)
    @test_throws ArgumentError estimate_crossings(
        MersenneTwister(1), scan; confidence=1.0)
end
