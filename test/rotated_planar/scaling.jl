function synthetic_channel_scan(; crossing::Union{Nothing,Float64}=0.1,
        nu::Float64=1.25, profile::Symbol=:identifiable,
        error_scale::Float64=1.0)
    distances = [3, 5, 7]
    error_rates = [0.04, 0.07, 0.1, 0.13, 0.16]
    shots = 100_000_000
    slopes = Dict(3 => 1.0, 5 => 0.8, 7 => 0.6)
    points = LogicalFailurePoint[]
    for distance_value in distances, error_rate in error_rates
        selected_rate = if crossing === nothing
            0.20 + 0.01 * ((distance_value - 1) ÷ 2) +
                0.5 * (error_rate - 0.1)
        elseif profile === :identifiable
            0.25 + 0.025 * (error_rate - crossing) * distance_value^(1 / nu)
        elseif profile === :boundary || profile === :flat
            0.25 + slopes[distance_value] * (error_rate - crossing)
        else
            error("unknown synthetic profile")
        end
        selected_count = round(Int, shots * selected_rate)
        selected_rate = selected_count / shots
        selected_error = profile === :flat ? 1e153 :
            error_scale * sqrt(selected_rate * (1 - selected_rate) / shots)
        push!(points, LogicalFailurePoint(
            distance_value, :x_ns, :as, :zero, :gate_layer,
            error_rate, 0.0, shots, 17,
            selected_count, selected_rate, selected_error,
            selected_count, selected_rate, selected_error,
            selected_count, selected_rate, selected_error,
            selected_count, selected_rate, selected_error,
        ))
    end
    return ChannelFailureScan(
        :as, :x_only, :logical_x, :zero, :x_ns, :gate_layer,
        distances, error_rates, shots, shots, 17, UInt64(17), points)
end

function with_logical_x_rate(point::LogicalFailurePoint, rate::Float64)
    count = round(Int, point.shots * rate)
    rate = count / point.shots
    standard_error = sqrt(rate * (1 - rate) / point.shots)
    return LogicalFailurePoint(
        point.distance, point.boundary_orientation, point.construction,
        point.logical_state, point.clock, point.p_x, point.p_z,
        point.shots, point.seed, count, rate, standard_error,
        point.logical_z_failures, point.logical_z_failure_rate,
        point.logical_z_standard_error, point.any_logical_failures,
        point.any_logical_failure_rate, point.any_logical_standard_error,
        point.state_failures, point.state_failure_rate,
        point.state_standard_error)
end

"""A scan with a zero-failure plateau and one spurious low-p crossing."""
function plateau_crossing_scan()
    distances = [3, 5, 7]
    error_rates = [0.0, 0.005, 0.01, 0.04, 0.06, 0.07, 0.1]
    rates = Dict(
        3 => [0.0, 0.0, 0.006, 0.019, 0.070, 0.125, 0.200],
        5 => [0.0, 0.0, 0.005, 0.020, 0.075, 0.120, 0.220],
        7 => [0.0, 0.0, 0.004, 0.015, 0.070, 0.125, 0.240],
    )
    shots = 1_000_000
    points = LogicalFailurePoint[]
    for distance_value in distances, (index, error_rate) in enumerate(error_rates)
        rate = rates[distance_value][index]
        count = round(Int, shots * rate)
        standard_error = sqrt(rate * (1 - rate) / shots)
        push!(points, LogicalFailurePoint(
            distance_value, :x_ns, :as, :zero, :gate_layer,
            error_rate, 0.0, shots, 23,
            count, rate, standard_error,
            count, rate, standard_error,
            count, rate, standard_error,
            count, rate, standard_error,
        ))
    end
    return ChannelFailureScan(
        :as, :x_only, :logical_x, :zero, :x_ns, :gate_layer,
        distances, error_rates, shots, shots, 23, UInt64(23), points)
end

@testset "Exploratory scaling fit" begin
    @test isdefined(TopoNoise, :PairCrossing)
    @test isdefined(TopoNoise, :ThresholdFit)
    @test isdefined(TopoNoise, :FittedConstructionChannelComparison)
    @test isdefined(TopoNoise, :fit_channel_threshold)
    @test isdefined(TopoNoise, :fit_construction_channel_comparison)
    @test isdefined(TopoNoise, :scaled_error_rate)

    state_point = LogicalFailurePoint(
        3, :x_ns, :as, :zero, :gate_layer, 0.0, 0.1, 10, 31,
        1, 0.1, sqrt(0.1 * 0.9 / 10),
        7, 0.7, sqrt(0.7 * 0.3 / 10),
        7, 0.7, sqrt(0.7 * 0.3 / 10),
        1, 0.1, sqrt(0.1 * 0.9 / 10),
    )
    state_resampled = TopoNoise._with_selected_channel_count(
        state_point, :state_failure, 2)
    @test state_resampled.logical_x_failures == 1
    @test state_resampled.logical_z_failures == 7
    @test state_resampled.state_failures == 2
    @test state_resampled.state_failure_rate == 0.2

    scan = synthetic_channel_scan(; crossing=0.1, nu=1.25)
    fit = fit_channel_threshold(scan; bootstrap_replicates=0, bootstrap_seed=7)
    @test fit.status === :success
    @test isapprox(fit.p_c, 0.1; atol=1e-10)
    @test isapprox(fit.nu, 1.25; atol=0.01)
    @test length(fit.crossings) == 2
    @test fit.bootstrap_replicates == 0
    @test fit.p_c_bootstrap_interval === nothing
    @test fit.nu_bootstrap_interval === nothing
    @test fit.exploratory
    @test scaled_error_rate(0.11, 5, fit) ==
        (0.11 - fit.p_c) * 5^(1 / fit.nu)

    endpoint_fit = fit_channel_threshold(
        synthetic_channel_scan(; crossing=0.1, profile=:boundary);
        bootstrap_replicates=0, bootstrap_seed=9)
    @test endpoint_fit.status === :unavailable
    @test occursin("boundary", endpoint_fit.diagnostic)

    flat_fit = fit_channel_threshold(
        synthetic_channel_scan(; crossing=0.1, profile=:flat);
        bootstrap_replicates=0, bootstrap_seed=9)
    @test flat_fit.status === :unavailable
    @test occursin("identifiable", flat_fit.diagnostic)

    unavailable = synthetic_channel_scan(; crossing=nothing)
    failed = fit_channel_threshold(unavailable; bootstrap_replicates=10, bootstrap_seed=8)
    @test failed.status === :unavailable
    @test failed.p_c === nothing
    @test failed.nu === nothing
    @test failed.bootstrap_successes == 0
    @test !isempty(failed.diagnostic)

    ambiguous = synthetic_channel_scan(; crossing=0.1)
    for index in eachindex(ambiguous.points)
        point = ambiguous.points[index]
        if point.distance == 3
            rate = [0.15, 0.30, 0.15, 0.30, 0.15][findfirst(==(point.p_x), ambiguous.error_rates)]
            ambiguous.points[index] = with_logical_x_rate(point, rate)
        elseif point.distance == 5
            ambiguous.points[index] = with_logical_x_rate(point, 0.20)
        end
    end
    ambiguous_fit = fit_channel_threshold(ambiguous; bootstrap_replicates=5)
    @test ambiguous_fit.status === :unavailable
    @test !isempty(ambiguous_fit.diagnostic)

    plateau_crossings, plateau_diagnostic = TopoNoise._adjacent_crossings(
        plateau_crossing_scan())
    @test plateau_diagnostic === nothing
    @test length(plateau_crossings) == 2
    @test all(crossing -> 0.06 < crossing.p < 0.07, plateau_crossings)

    first = fit_channel_threshold(scan; bootstrap_replicates=20, bootstrap_seed=42)
    second = fit_channel_threshold(scan; bootstrap_replicates=20, bootstrap_seed=42)
    @test (
        first.status, first.p_c, first.p_c_standard_error, first.nu,
        first.nu_standard_error, first.bootstrap_successes,
        first.p_c_bootstrap_interval, first.nu_bootstrap_interval,
    ) == (
        second.status, second.p_c, second.p_c_standard_error, second.nu,
        second.nu_standard_error, second.bootstrap_successes,
        second.p_c_bootstrap_interval, second.nu_bootstrap_interval,
    )
    @test first.status === :success
    @test first.bootstrap_successes >= 2
    @test first.p_c_bootstrap_interval !== nothing
    @test first.nu_bootstrap_interval !== nothing

    insufficient_bootstrap = fit_channel_threshold(
        scan; bootstrap_replicates=1, bootstrap_seed=42)
    @test insufficient_bootstrap.status === :unavailable
    @test occursin("bootstrap", insufficient_bootstrap.diagnostic)

    synthetic_series = [
        ChannelFailureScan(
            construction, channel,
            channel === :x_only ? :logical_x : :logical_z, :zero, :x_ns,
            :gate_layer, copy(scan.distances), copy(scan.error_rates), scan.shots,
            scan.batch_size, 3, UInt64(index), copy(scan.points),
        )
        for (index, (construction, channel)) in
            enumerate(((:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)))
    ]
    raw = ConstructionChannelComparison(
        :zero, :x_ns, :gate_layer, copy(scan.distances), copy(scan.error_rates),
        scan.shots, scan.batch_size, 3, synthetic_series,
    )
    fitted = fit_construction_channel_comparison(
        MersenneTwister(3), raw; bootstrap_replicates=10)
    @test length(fitted.fits) == 4
    @test [(fit.construction, fit.error_channel) for fit in fitted.fits] ==
        [(:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)]
    @test all(fit -> fit.status === :success, fitted.fits)
end
