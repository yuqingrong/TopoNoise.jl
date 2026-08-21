using Random
using Test

function _planar_zero_errors(model)
    IsometricPlanarVirtualErrors(
        falses(model.rows, model.cols - 1), falses(model.rows - 1, model.cols))
end

function _single_planar_error(model, direction::Symbol, row::Int, col::Int)
    errors = _planar_zero_errors(model)
    getproperty(errors, direction)[row, col] = true
    return errors
end

@testset "Pseudo-toric API is retired" begin
    @test !isdefined(TopoNoise, :IsometricToricCodeModel)
    @test !isdefined(TopoNoise, :scan_isometric_toric_capacity)
end

@testset "Isometric planar-code geometry and detectors" begin
    model = IsometricPlanarCodeModel(3)
    @test model.distance == 3
    @test (model.rows, model.cols) == (3, 3)

    errors = IsometricPlanarVirtualErrors(falses(3, 2), falses(2, 3))
    syndrome = isometric_planar_syndrome(errors)
    @test size(syndrome.plaquettes) == (2, 2)
    @test length(syndrome.north) == 2
    @test length(syndrome.south) == 2
    @test !any(syndrome.plaquettes)
    @test !any(syndrome.north)
    @test !any(syndrome.south)
end

@testset "Isometric planar-code Yao detector reference" begin
    for distance in (2, 3)
        model = IsometricPlanarCodeModel(distance)
        error_records = [_planar_zero_errors(model)]
        for row in 1:model.rows, col in 1:(model.cols - 1)
            push!(error_records, _single_planar_error(model, :horizontal, row, col))
        end
        for row in 1:(model.rows - 1), col in 1:model.cols
            push!(error_records, _single_planar_error(model, :vertical, row, col))
        end
        append!(error_records, [
            sample_isometric_planar_virtual_errors(
                MersenneTwister(1000 + distance + sample), model; error_rate=0.2)
            for sample in 1:3])

        for (index, errors) in enumerate(error_records)
            trajectory = sample_isometric_planar_yao_trajectory(
                MersenneTwister(2000 + 100 * distance + index), model, errors)
            yao_syndrome = isometric_planar_syndrome(
                isometric_planar_bond_mismatches(trajectory))
            direct_syndrome = isometric_planar_syndrome(errors)
            @test yao_syndrome.plaquettes == direct_syndrome.plaquettes
            @test yao_syndrome.north == direct_syndrome.north
            @test yao_syndrome.south == direct_syndrome.south
            @test TopoNoise._isometric_planar_yao_logical_frame(trajectory) ==
                  TopoNoise._isometric_planar_logical_frame(errors)
        end
    end
end

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

@testset "Isometric planar-code encoder and parity-check layout" begin
    model = IsometricPlanarCodeModel(3)
    encoder = isometric_planar_encoder(model)
    @test encoder.wire_count == 36
    @test size(encoder.output_wires) == (3, 3, 4)
    @test sort(vec(encoder.output_wires)) == collect(1:36)

    checks = isometric_planar_check_layer(encoder)
    @test length(checks.detector_supports) == 8
    @test all(length(support) == 8 for support in checks.detector_supports[1:4])
    @test all(length(support) == 2 for support in checks.detector_supports[5:8])
    @test length(checks.logical_z_support) == 4
end

@testset "Isometric planar-code topology and decoder" begin
    for distance in 2:8
        model = IsometricPlanarCodeModel(distance)
        @test TopoNoise._isometric_planar_minimum_odd_logical_weight(model) == distance
    end

    model = IsometricPlanarCodeModel(4)
    logical = _planar_zero_errors(model)
    logical.vertical[1, :] .= true
    logical_syndrome = isometric_planar_syndrome(logical)
    @test !any(logical_syndrome.plaquettes)
    @test !any(logical_syndrome.north)
    @test !any(logical_syndrome.south)
    @test isometric_planar_logical_failure(logical, Correction(4, 4))

    corner = _planar_zero_errors(model)
    corner.horizontal[1, 1] = true
    corner.vertical[1, 1] = true
    corner_syndrome = isometric_planar_syndrome(corner)
    @test corner_syndrome.north == BitVector([true, false, false])
    @test !any(corner_syndrome.plaquettes)

    for distance in (3, 4)
        model = IsometricPlanarCodeModel(distance)
        for row in 1:model.rows, col in 1:(model.cols - 1)
            errors = _single_planar_error(model, :horizontal, row, col)
            correction = decode_isometric_planar_syndrome(
                model, isometric_planar_syndrome(errors))
            @test !isometric_planar_logical_failure(errors, correction)
        end
        for row in 1:(model.rows - 1), col in 1:model.cols
            errors = _single_planar_error(model, :vertical, row, col)
            correction = decode_isometric_planar_syndrome(
                model, isometric_planar_syndrome(errors))
            @test !isometric_planar_logical_failure(errors, correction)
        end
    end

    diagnostic = IsometricPlanarCodeModel(2)
    west = _single_planar_error(diagnostic, :vertical, 1, 1)
    correction = decode_isometric_planar_syndrome(
        diagnostic, isometric_planar_syndrome(west))
    @test correction.vertical[1, 1]
end

@testset "Isometric planar-code distance-five low-weight decoding" begin
    model = IsometricPlanarCodeModel(5)
    supports = Tuple{Symbol,Int,Int}[]
    append!(supports, [(:horizontal, row, col)
                       for row in 1:model.rows, col in 1:(model.cols - 1)])
    append!(supports, [(:vertical, row, col)
                       for row in 1:(model.rows - 1), col in 1:model.cols])

    records = IsometricPlanarVirtualErrors[]
    for support in supports
        errors = _planar_zero_errors(model)
        getproperty(errors, support[1])[support[2], support[3]] = true
        push!(records, errors)
    end
    for first in 1:(length(supports) - 1), second in (first + 1):length(supports)
        errors = _planar_zero_errors(model)
        for support in (supports[first], supports[second])
            getproperty(errors, support[1])[support[2], support[3]] = true
        end
        push!(records, errors)
    end

    detector_batch = falses(length(records), 24)
    scalar_frames = falses(length(records))
    for (index, errors) in enumerate(records)
        syndrome = isometric_planar_syndrome(errors)
        detector_batch[index, :] .= TopoNoise._isometric_planar_syndrome_vector(
            model, syndrome) .!= 0
        correction = decode_isometric_planar_syndrome(model, syndrome)
        @test !isometric_planar_logical_failure(errors, correction)
        scalar_frames[index] = xor(correction.vertical[:, 1]...)
    end
    @test TopoNoise._isometric_planar_batch_predictions(model, detector_batch) == scalar_frames
end
