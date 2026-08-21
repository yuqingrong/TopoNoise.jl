using Random
using Test

@testset "Isometric planar-code batched capacity estimates" begin
    model = IsometricPlanarCodeModel(3)
    zero_point = estimate_isometric_planar_capacity(
        MersenneTwister(10), model, 0.0; shots=20, batches=2)
    @test zero_point.distance == 3
    @test zero_point.logical_failure_count == 0
    @test zero_point.logical_failure_rate == 0.0
    @test length(zero_point.logical_failure_batches) == 2

    point = estimate_isometric_planar_capacity(
        MersenneTwister(11), model, 0.1; shots=20, batches=2)
    @test 0.0 <= point.logical_failure_rate <= 1.0
    @test length(point.logical_failure_batches) == 2

    scan = scan_isometric_planar_capacity(
        MersenneTwister(12), [3, 5], [0.0, 0.1]; shots=20, batches=2)
    @test scan.distances == [3, 5]
    @test length(scan.points) == 4
    @test length(estimate_isometric_planar_crossings(
        MersenneTwister(13), scan; bootstrap=8)) == 1

    rates = [0.05, 0.10, 0.15, 0.20]
    @test TopoNoise._unique_crossing(
        rates, [0.1, 0.3, 0.7, 0.9], [0.2, 0.4, 0.6, 0.8]).status == :ok
    @test TopoNoise._unique_crossing(
        rates, [0.1, 0.8, 0.2, 0.9], [0.2, 0.3, 0.7, 0.4]).status == :unstable
end

function _synthetic_isometric_planar_scaling_scan(distances=[5, 7, 9, 11])
    rates = collect(0.07:0.01:0.13)
    points = IsometricPlanarCapacityPoint[]
    for distance in distances, rate in rates
        x = (rate - 0.10) * distance^(1 / 1.5)
        failure_rate = 0.25 + 0.60 * x + 0.20 * x^3
        batches = [failure_rate - 0.004, failure_rate - 0.002,
                   failure_rate + 0.002, failure_rate + 0.004]
        push!(points, IsometricPlanarCapacityPoint(
            distance, rate, 400, round(Int, 400 * failure_rate), failure_rate,
            0.002, batches, distance == 2, 17))
    end
    return IsometricPlanarCapacityScan(distances, rates, points)
end

function _diagnostic_isometric_planar_scaling_scan(;
        distances=[5, 7, 9, 11, 13], p_c=0.1037, nu=1.43,
        flat=false, zero_batches=false, failed_batches=false)
    rates = collect(0.07:0.01:0.14)
    offsets = [-0.012, -0.009, -0.005, -0.002, 0.003, 0.006, 0.008, 0.011]
    points = IsometricPlanarCapacityPoint[]
    for (distance_index, distance) in enumerate(distances),
            (rate_index, rate) in enumerate(rates)
        x = (rate - p_c) * distance^(1 / nu)
        failure_rate = flat ? 0.25 : 0.25 + 0.72x + 0.18x^2 - 0.10x^3
        batch_scale = 0.65 + 0.04 * distance_index + 0.02 * rate_index
        batches = if failed_batches
            fill(NaN, length(offsets))
        elseif zero_batches
            fill(failure_rate, length(offsets))
        else
            failure_rate .+ batch_scale .* circshift(offsets, rate_index)
        end
        push!(points, IsometricPlanarCapacityPoint(
            distance, rate, 800, round(Int, 800 * failure_rate), failure_rate,
            0.002, batches, false, 29))
    end
    return IsometricPlanarCapacityScan(collect(distances), rates, points)
end

@testset "Isometric planar finite-size scaling fit" begin
    @test isdefined(TopoNoise, :IsometricPlanarScalingFit)
    @test isdefined(TopoNoise, :fit_isometric_planar_scaling)

    if isdefined(TopoNoise, :fit_isometric_planar_scaling)
        scan = _synthetic_isometric_planar_scaling_scan()
        fit = fit_isometric_planar_scaling(
            MersenneTwister(23), scan; bootstrap=24)
        @test fit.status == :ok
        @test isapprox(fit.p_c, 0.10; atol=0.005)
        @test isapprox(fit.nu, 1.5; atol=0.15)
        @test fit.p_c_se >= 0
        @test fit.nu_se >= 0
        @test fit.valid_bootstrap_fraction >= 0.8

        diagnostic_scan = _synthetic_isometric_planar_scaling_scan([2, 5, 7, 9, 11])
        diagnostic_fit = fit_isometric_planar_scaling(
            MersenneTwister(24), diagnostic_scan; bootstrap=24)
        @test diagnostic_fit.status == :ok
        @test isapprox(diagnostic_fit.p_c, fit.p_c; atol=0.005)

        short_scan = _synthetic_isometric_planar_scaling_scan([5, 7])
        short_fit = fit_isometric_planar_scaling(
            MersenneTwister(25), short_scan; bootstrap=24)
        @test short_fit.status == :insufficient_sizes
    end
end


@testset "Isometric planar scaling diagnostics" begin
    @test isdefined(TopoNoise, :IsometricPlanarScalingDiagnostics)
    @test isdefined(TopoNoise, :diagnose_isometric_planar_scaling)

    if isdefined(TopoNoise, :diagnose_isometric_planar_scaling)
        scan = _diagnostic_isometric_planar_scaling_scan()
        diagnostics = diagnose_isometric_planar_scaling(
            MersenneTwister(31), scan; bootstrap=32,
            loss_grid_size=(19, 17), sensitivity=true)

        @test diagnostics.fit.status == :ok
        @test isapprox(diagnostics.fit.p_c, 0.1037; atol=0.002)
        @test isapprox(diagnostics.fit.nu, 1.43; atol=0.12)
        @test size(diagnostics.loss_values) == (19, 17)
        @test all(isfinite, diagnostics.loss_values)
        minimum_index = argmin(diagnostics.loss_values)
        @test abs(diagnostics.loss_p_c[minimum_index[1]] - diagnostics.fit.p_c) <=
              step(range(first(diagnostics.loss_p_c), last(diagnostics.loss_p_c);
                         length=length(diagnostics.loss_p_c)))
        @test abs(diagnostics.loss_nu[minimum_index[2]] - diagnostics.fit.nu) <=
              step(range(first(diagnostics.loss_nu), last(diagnostics.loss_nu);
                         length=length(diagnostics.loss_nu)))
        @test length(diagnostics.bootstrap_samples) == 32
        successful = filter(sample -> sample.status == :ok,
                            diagnostics.bootstrap_samples)
        @test length(successful) >= 26
        @test maximum(sample.p_c for sample in successful) >
              minimum(sample.p_c for sample in successful)
        @test maximum(sample.nu for sample in successful) >
              minimum(sample.nu for sample in successful)
        @test any(sample.input_rms_shift > 0 for sample in successful)
        @test all(sample.evaluations > 0 for sample in successful)

        rates = [point.error_rate for point in
                 TopoNoise._isometric_planar_scaling_selected_points(
                     scan, (0.05, 0.45))]
        p_grid = collect(range(minimum(rates), maximum(rates); length=13))
        @test minimum(abs(diagnostics.fit.p_c - candidate) for candidate in p_grid) > 1e-5
        @test !isempty(diagnostics.optimizer_p_c)
        @test diagnostics.nominal_converged
        @test all(isfinite, diagnostics.optimizer_loss)
        @test diagnostics.optimizer_p_c[end] == diagnostics.fit.p_c
        @test diagnostics.optimizer_nu[end] == diagnostics.fit.nu
        @test any(summary.std_ratio > 0 for summary in diagnostics.batch_diagnostics)
        @test any(summary.unique_batch_rates > 1
                  for summary in diagnostics.batch_diagnostics)
        @test Set(summary.label for summary in diagnostics.sensitivities) == Set([
            "baseline", "drop d=5", "drop d=7", "drop d=9", "drop d=11",
            "drop d=13", "window 0.05-0.35", "window 0.10-0.45"])

        flat = diagnose_isometric_planar_scaling(
            MersenneTwister(32),
            _diagnostic_isometric_planar_scaling_scan(flat=true);
            bootstrap=8, loss_grid_size=(9, 7), sensitivity=false)
        @test maximum(flat.loss_values) - minimum(flat.loss_values) < 1e-8

        zero_batches = diagnose_isometric_planar_scaling(
            MersenneTwister(33),
            _diagnostic_isometric_planar_scaling_scan(zero_batches=true);
            bootstrap=8, loss_grid_size=(9, 7), sensitivity=false)
        @test all(summary.std_ratio <= 1e-12
                  for summary in zero_batches.batch_diagnostics)
        @test all(summary.unique_batch_rates == 1
                  for summary in zero_batches.batch_diagnostics)
        @test all(sample.input_rms_shift <= 1e-12
                  for sample in zero_batches.bootstrap_samples)

        failed = diagnose_isometric_planar_scaling(
            MersenneTwister(34),
            _diagnostic_isometric_planar_scaling_scan(failed_batches=true);
            bootstrap=8, loss_grid_size=(9, 7), sensitivity=false)
        @test failed.fit.status == :unstable
        @test all(sample.status == :fit_failed for sample in failed.bootstrap_samples)

        boundary = diagnose_isometric_planar_scaling(
            MersenneTwister(35),
            _diagnostic_isometric_planar_scaling_scan(p_c=0.05);
            bootstrap=8, loss_grid_size=(9, 7), sensitivity=false)
        @test boundary.nominal_boundary_hit

        unavailable = diagnose_isometric_planar_scaling(
            MersenneTwister(36),
            _diagnostic_isometric_planar_scaling_scan(distances=[5, 7]);
            bootstrap=8, loss_grid_size=(9, 7), sensitivity=true)
        @test unavailable.fit.status == :insufficient_sizes
        @test isempty(unavailable.bootstrap_samples)
        @test isempty(unavailable.optimizer_p_c)
    end
end
