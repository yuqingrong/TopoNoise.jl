@testset "Rotated-planar construction channel comparisons" begin
    comparison_api = (
        :ChannelFailureScan,
        :ConstructionChannelComparison,
        :scan_channel_logical_failure,
        :run_construction_channel_comparison,
        :comparison_series,
        :channel_failure_count,
        :channel_failure_rate,
        :channel_failure_standard_error,
    )
    for name in comparison_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), comparison_api)
        @testset "zero-noise channel scans select their intended observable" begin
            x_scan = scan_channel_logical_failure(
                MersenneTwister(21);
                distances=[3], error_rates=[0.0], construction=:as,
                error_channel=:x_only, shots=2, batch_size=2, seed=21)
            z_scan = scan_channel_logical_failure(
                MersenneTwister(22);
                distances=[3], error_rates=[0.0], construction=:bp,
                error_channel=:z_only, shots=2, batch_size=2, seed=22)

            @test x_scan isa ChannelFailureScan
            @test x_scan.logical_state === :plus
            @test x_scan.logical_observable === :logical_x
            @test z_scan.logical_state === :zero
            @test z_scan.logical_observable === :logical_z
            for scan in (x_scan, z_scan)
                point = only(scan.points)
                @test point.p_x == 0.0
                @test point.p_z == 0.0
                @test channel_failure_count(scan, point) == 0
                @test channel_failure_rate(scan, point) == 0.0
                @test channel_failure_standard_error(scan, point) == 0.0
            end

            @test_throws ArgumentError scan_channel_logical_failure(
                MersenneTwister(23);
                distances=[3], error_rates=[0.0], error_channel=:balanced,
                shots=2, batch_size=2)
        end

        @testset "state-failure scans select the prepared-state metric" begin
            state_scan = scan_channel_logical_failure(
                MersenneTwister(24);
                distances=[3], error_rates=[0.0], construction=:bp,
                error_channel=:z_only, failure_metric=:state_failure,
                shots=2, batch_size=2, seed=24)
            @test state_scan.logical_observable === :state_failure

            point = LogicalFailurePoint(
                3, :x_ns, :bp, :zero, :plaquette, 0.0, 0.1, 10, 24,
                1, 0.1, sqrt(0.1 * 0.9 / 10),
                7, 0.7, sqrt(0.7 * 0.3 / 10),
                7, 0.7, sqrt(0.7 * 0.3 / 10),
                1, 0.1, sqrt(0.1 * 0.9 / 10),
            )
            selector = ChannelFailureScan(
                :bp, :z_only, :state_failure, :zero, :x_ns, :plaquette,
                [3], [0.1], 10, 10, 24, UInt64(24), [point])
            @test channel_failure_count(selector, point) == 1
            @test channel_failure_rate(selector, point) == 0.1
            @test channel_failure_standard_error(selector, point) ==
                  sqrt(0.1 * 0.9 / 10)
            @test_throws ArgumentError scan_channel_logical_failure(
                MersenneTwister(25);
                distances=[3], error_rates=[0.0], error_channel=:z_only,
                failure_metric=:invalid, shots=2, batch_size=2)
        end

        @testset "comparison ordering and independently seeded series are reproducible" begin
            comparison = run_construction_channel_comparison(
                MersenneTwister(99);
                distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=99)
            @test comparison isa ConstructionChannelComparison
            @test comparison.logical_state === :native
            @test [scan.logical_state for scan in comparison.series] ==
                  [:plus, :plus, :zero, :zero]
            @test all(point.logical_state === scan.logical_state
                      for scan in comparison.series for point in scan.points)
            @test [(scan.construction, scan.error_channel) for scan in comparison.series] == [
                (:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only),
            ]
            @test comparison_series(comparison, :bp, :z_only) === comparison.series[4]
            @test_throws ArgumentError comparison_series(comparison, :bp, :balanced)

            x_scan = scan_channel_logical_failure(
                MersenneTwister(41);
                distances=[3], error_rates=[0.1], construction=:as,
                error_channel=:x_only, shots=4, batch_size=2, seed=41)
            z_scan = scan_channel_logical_failure(
                MersenneTwister(42);
                distances=[3], error_rates=[0.1], construction=:bp,
                error_channel=:z_only, shots=4, batch_size=2, seed=42)
            @test only(x_scan.points).p_x == 0.1
            @test only(x_scan.points).p_z == 0.0
            @test only(z_scan.points).p_x == 0.0
            @test only(z_scan.points).p_z == 0.1

            first = run_construction_channel_comparison(
                MersenneTwister(99);
                distances=[3], error_rates=[0.1], shots=4, batch_size=2, seed=99)
            repeated = run_construction_channel_comparison(
                MersenneTwister(99);
                distances=[3], error_rates=[0.1], shots=4, batch_size=2, seed=99)
            @test [scan.series_seed for scan in first.series] ==
                  [scan.series_seed for scan in repeated.series]
            project(scan) = [
                (scan.construction, scan.error_channel, scan.logical_observable,
                 point.p_x, point.p_z, point.logical_x_failures,
                 point.logical_z_failures)
                for point in scan.points
            ]
            @test [project(scan) for scan in first.series] ==
                  [project(scan) for scan in repeated.series]
        end

        @testset "native preparation reports logical errors even when the state is unchanged" begin
            for (construction, channel, state, observable) in
                ((:as, :x_only, :plus, :logical_x), (:bp, :z_only, :zero, :logical_z))
                scan = scan_channel_logical_failure(MersenneTwister(72);
                    distances=[3], error_rates=[0.15], construction,
                    error_channel=channel, shots=512, batch_size=256)
                point = only(scan.points)
                @test scan.logical_state === point.logical_state === state
                @test scan.logical_observable === observable
                @test point.state_failures == 0
                @test channel_failure_count(scan, point) > 0
                @test channel_failure_rate(scan, point) ==
                      channel_failure_count(scan, point) / point.shots
            end
        end

        @testset "plus-state Z scans measure phase-sensitive state failure" begin
            scan = scan_channel_logical_failure(
                MersenneTwister(73);
                distances=[3], error_rates=[0.0, 0.12], construction=:as,
                logical_state=:plus, error_channel=:z_only,
                failure_metric=:state_failure, shots=128, batch_size=64)
            @test scan.logical_state === :plus
            @test all(point.logical_state === :plus for point in scan.points)
            @test first(scan.points).state_failures == 0
            noisy = last(scan.points)
            @test noisy.logical_x_failures == 0
            @test noisy.logical_z_failures > 0
            @test noisy.state_failures == noisy.logical_z_failures
            @test channel_failure_count(scan, noisy) == noisy.logical_z_failures

            comparison = run_construction_channel_comparison(
                MersenneTwister(74); logical_state=:plus,
                distances=[3], error_rates=[0.0], shots=2, batch_size=2)
            @test comparison.logical_state === :plus
            @test all(series.logical_state === :plus for series in comparison.series)
            @test all(only(series.points).logical_state === :plus
                      for series in comparison.series)
            @test_throws ArgumentError scan_channel_logical_failure(
                MersenneTwister(75); logical_state=:invalid,
                distances=[3], error_rates=[0.0], error_channel=:z_only, shots=2)
        end

        @testset "keyed series seeds are independent of collection order" begin
            root_seed = UInt64(0x123456789abcdef0)
            canonical = Dict(
                key => TopoNoise._comparison_series_seed(root_seed, key)
                for key in ((:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only))
            )
            reordered = Dict(
                key => TopoNoise._comparison_series_seed(root_seed, key)
                for key in reverse(((:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)))
            )
            extended = Dict(
                key => TopoNoise._comparison_series_seed(root_seed, key)
                for key in ((:unused, :other), (:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only))
            )
            @test canonical == reordered
            @test all(key -> canonical[key] == extended[key], keys(canonical))
            @test length(unique(values(canonical))) == 4
        end
    end
end
