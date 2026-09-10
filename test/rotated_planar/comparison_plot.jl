using CairoMakie

function synthetic_fitted_comparison()
    distances = [3, 5, 7]
    error_rates = [0.05, 0.10, 0.15]
    shots = 10_000
    series = ChannelFailureScan[]
    fits = ThresholdFit[]
    combinations = ((:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only))

    for (index, (construction, error_channel)) in enumerate(combinations)
        observable = error_channel === :x_only ? :logical_x : :logical_z
        points = LogicalFailurePoint[]
        for distance_value in distances, error_rate in error_rates
            rate = 0.25 + (8 - distance_value) / 10 * (error_rate - 0.10)
            count = round(Int, shots * rate)
            rate = count / shots
            standard_error = sqrt(rate * (1 - rate) / shots)
            p_x, p_z = error_channel === :x_only ? (error_rate, 0.0) : (0.0, error_rate)
            push!(points, LogicalFailurePoint(
                distance_value, :x_ns, construction, :zero, :gate_layer,
                p_x, p_z, shots, 1234,
                count, rate, standard_error,
                count, rate, standard_error,
                count, rate, standard_error,
                count, rate, standard_error,
            ))
        end
        scan = ChannelFailureScan(
            construction, error_channel, observable, :zero, :x_ns, :gate_layer,
            copy(distances), copy(error_rates), shots, shots, 1234, UInt64(index), points)
        push!(series, scan)
        crossings = PairCrossing[
            PairCrossing(3, 5, 0.10, 0.002),
            PairCrossing(5, 7, 0.10, 0.003),
        ]
        push!(fits, ThresholdFit(
            construction, error_channel, observable, :success,
            "exploratory \"synthetic\" fit", crossings, 0.10, 0.002,
            1.25, 0.10, 20, 19, UInt64(42 + index), true))
    end
    raw = ConstructionChannelComparison(
        :zero, :x_ns, :gate_layer, distances, error_rates, shots, shots, 1234, series)
    return FittedConstructionChannelComparison(raw, fits)
end

@testset "Construction/channel comparison artifacts" begin
    artifact_api = (
        :write_comparison_raw_csv,
        :write_comparison_fit_csv,
        :plot_construction_channel_comparison,
        :plot_construction_channel_panel,
        :save_construction_channel_comparison,
    )
    for name in artifact_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), artifact_api)
        fitted = synthetic_fitted_comparison()
        mktempdir() do directory
            paths = save_construction_channel_comparison(fitted, directory; basename="comparison")
            @test keys(paths) == (:raw_csv, :fits_csv, :combined, :panels)
            @test isfile(paths.raw_csv)
            @test isfile(paths.fits_csv)
            @test readlines(paths.raw_csv)[1] ==
                "construction,error_channel,logical_observable,distance,p_x,p_z,shots,master_seed,series_seed,logical_failures,logical_failure_rate,logical_failure_standard_error,logical_state,boundary_orientation,clock"
            @test readlines(paths.fits_csv)[1] ==
                "construction,error_channel,logical_observable,status,diagnostic,crossing_lower_distance,crossing_upper_distance,crossing_p,crossing_standard_error,p_c,p_c_standard_error,nu,nu_standard_error,bootstrap_replicates,bootstrap_successes,bootstrap_seed,exploratory"
            @test length(readlines(paths.raw_csv)) == 1 + 4 * 3 * 3
            @test length(readlines(paths.fits_csv)) == 1 + 4 * 2
            @test length(paths.panels) == 4
            for group in (paths.combined, values(paths.panels)...), path in values(group)
                @test isfile(path)
                @test filesize(path) > 100
            end
            combined = plot_construction_channel_comparison(fitted)
            @test combined isa CairoMakie.Figure
            @test count(
                item -> item isa CairoMakie.Axis && item.xlabel[] == "Physical error rate",
                combined.content) == 4
            @test plot_construction_channel_panel(fitted, :as, :x_only) isa CairoMakie.Figure
            @test_throws ArgumentError plot_construction_channel_panel(fitted, :invalid, :x_only)
        end

        unavailable = deepcopy(fitted)
        unavailable.fits[1] = ThresholdFit(
            :as, :x_only, :logical_x, :unavailable,
            "missing \"crossing\"", PairCrossing[], nothing, nothing,
            nothing, nothing, 0, 0, UInt64(0), true)
        mktempdir() do directory
            raw_path = joinpath(directory, "raw.csv")
            fits_path = joinpath(directory, "fits.csv")
            write_comparison_raw_csv(raw_path, unavailable)
            write_comparison_fit_csv(fits_path, unavailable)
            raw_rows = readlines(raw_path)
            fit_rows = readlines(fits_path)
            @test any(
                row -> occursin(",0.05,0.0,", row),
                filter(row -> startswith(row, "as,x_only"), raw_rows[2:end]))
            unavailable_row = only(filter(row -> startswith(row, "as,x_only,logical_x,unavailable"), fit_rows[2:end]))
            @test occursin("\"missing \"\"crossing\"\"\"", unavailable_row)
            @test occursin(",,,,,", unavailable_row)
            @test plot_construction_channel_panel(unavailable, :as, :x_only) isa CairoMakie.Figure
        end

        mktempdir() do directory
            for basename in ("", ".", "..", "../escape", "a/b", "a\\b")
                @test_throws ArgumentError save_construction_channel_comparison(
                    fitted, directory; basename)
            end
        end
    end
end

const _COMPARISON_CLI = joinpath(
    @__DIR__, "..", "..", "examples", "compare_rotated_planar_constructions.jl")

@testset "Construction/channel comparison CLI" begin
    @test isfile(_COMPARISON_CLI)
    if isfile(_COMPARISON_CLI)
        include(_COMPARISON_CLI)
        mktempdir() do directory
            output, errors = IOBuffer(), IOBuffer()
            result = RotatedPlanarConstructionComparison.main([
                "--distances", "3", "--error-rates", "0", "--shots", "2",
                "--batch-size", "2", "--bootstrap-replicates", "2", "--no-fit",
                "--seed", "5", "--output-dir", directory, "--basename", "smoke",
            ]; io=output, error_io=errors)
            @test result == 0
            @test isfile(joinpath(directory, "smoke-combined.svg"))
            @test isfile(joinpath(directory, "smoke-as-x-only.png"))
            @test contains(String(take!(output)), "fit unavailable")
            @test isempty(String(take!(errors)))
        end
        for arguments in (
                ["--bootstrap-replicates", "0"],
                ["--error-rates", "0.1,0.1"],
                ["--unknown"],
                ["--basename", "../escape"],
            )
            errors = IOBuffer()
            @test RotatedPlanarConstructionComparison.main(arguments; error_io=errors) == 1
            @test contains(String(take!(errors)), "error:")
        end
    end
end
