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
            @test x_scan.logical_state === :zero
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

        @testset "comparison ordering and independently seeded series are reproducible" begin
            comparison = run_construction_channel_comparison(
                MersenneTwister(99);
                distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=99)
            @test comparison isa ConstructionChannelComparison
            @test comparison.logical_state === :zero
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
    end
end
