const _COMPARISON_RAW_CSV_COLUMNS = (
    :construction,
    :error_channel,
    :logical_observable,
    :distance,
    :p_x,
    :p_z,
    :shots,
    :master_seed,
    :series_seed,
    :logical_failures,
    :logical_failure_rate,
    :logical_failure_standard_error,
    :logical_state,
    :boundary_orientation,
    :clock,
)

const _COMPARISON_FIT_CSV_COLUMNS = (
    :construction,
    :error_channel,
    :logical_observable,
    :status,
    :diagnostic,
    :crossing_lower_distance,
    :crossing_upper_distance,
    :crossing_p,
    :crossing_standard_error,
    :p_c,
    :p_c_standard_error,
    :p_c_bootstrap_interval_lower,
    :p_c_bootstrap_interval_upper,
    :nu,
    :nu_standard_error,
    :nu_bootstrap_interval_lower,
    :nu_bootstrap_interval_upper,
    :bootstrap_replicates,
    :bootstrap_successes,
    :bootstrap_seed,
    :exploratory,
)

_comparison_csv_value(value) = value === nothing ? "" : string(value)
_comparison_csv_diagnostic(value::AbstractString) =
    "\"" * replace(value, "\"" => "\"\"") * "\""

"""Write every selected-channel Monte Carlo point as a tidy CSV file."""
function write_comparison_raw_csv(
        path::AbstractString, fitted::FittedConstructionChannelComparison)
    output_path = abspath(path)
    open(output_path, "w") do io
        println(io, join(string.(_COMPARISON_RAW_CSV_COLUMNS), ','))
        for scan in fitted.raw.series, point in scan.points
            row = (
                scan.construction,
                scan.error_channel,
                scan.logical_observable,
                point.distance,
                point.p_x,
                point.p_z,
                point.shots,
                scan.master_seed,
                scan.series_seed,
                channel_failure_count(scan, point),
                channel_failure_rate(scan, point),
                channel_failure_standard_error(scan, point),
                scan.logical_state,
                scan.boundary_orientation,
                scan.clock,
            )
            println(io, join(_comparison_csv_value.(row), ','))
        end
    end
    return output_path
end

function _fit_csv_row(fit::ThresholdFit, crossing::Union{Nothing,PairCrossing})
    p_c_interval = fit.p_c_bootstrap_interval
    nu_interval = fit.nu_bootstrap_interval
    return (
        fit.construction,
        fit.error_channel,
        fit.logical_observable,
        fit.status,
        _comparison_csv_diagnostic(fit.diagnostic),
        crossing === nothing ? nothing : crossing.lower_distance,
        crossing === nothing ? nothing : crossing.upper_distance,
        crossing === nothing ? nothing : crossing.p,
        crossing === nothing ? nothing : crossing.standard_error,
        fit.p_c,
        fit.p_c_standard_error,
        p_c_interval === nothing ? nothing : first(p_c_interval),
        p_c_interval === nothing ? nothing : last(p_c_interval),
        fit.nu,
        fit.nu_standard_error,
        nu_interval === nothing ? nothing : first(nu_interval),
        nu_interval === nothing ? nothing : last(nu_interval),
        fit.bootstrap_replicates,
        fit.bootstrap_successes,
        fit.bootstrap_seed,
        fit.exploratory,
    )
end

"""Write panel fit summaries and their pairwise-crossing evidence as CSV."""
function write_comparison_fit_csv(
        path::AbstractString, fitted::FittedConstructionChannelComparison)
    output_path = abspath(path)
    open(output_path, "w") do io
        println(io, join(string.(_COMPARISON_FIT_CSV_COLUMNS), ','))
        for fit in fitted.fits
            crossings = isempty(fit.crossings) ? Union{Nothing,PairCrossing}[nothing] : fit.crossings
            for crossing in crossings
                row = _fit_csv_row(fit, crossing)
                values = map(enumerate(row)) do (index, value)
                    _COMPARISON_FIT_CSV_COLUMNS[index] === :diagnostic ? value :
                        _comparison_csv_value(value)
                end
                println(io, join(values, ','))
            end
        end
    end
    return output_path
end

function _comparison_fit(
        fitted::FittedConstructionChannelComparison,
        construction::Symbol,
        error_channel::Symbol)
    candidates = filter(
        fit -> fit.construction === construction && fit.error_channel === error_channel,
        fitted.fits)
    length(candidates) == 1 || throw(ArgumentError(
        "fitted comparison must contain exactly one fit for " *
        "construction=$(construction), error_channel=$(error_channel)"))
    return only(candidates)
end

function _channel_label(channel::Symbol)
    channel === :x_only && return "X-only noise"
    channel === :z_only && return "Z-only noise"
    throw(ArgumentError("unknown error channel: $(channel)"))
end

function _construction_label(construction::Symbol)
    construction === :as && return "As construction"
    construction === :bp && return "Bp construction"
    throw(ArgumentError("unknown construction: $(construction)"))
end

_physical_error_rate(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.error_channel === :x_only ? point.p_x : point.p_z

function _comparison_state_label(state::Symbol)
    label = state === :zero ? "0" : state === :one ? "1" :
            state === :plus ? "+" : state === :minus ? "-" : string(state)
    return "|$(label)_L>"
end

function _logical_axis_label(observable::Symbol, logical_state::Symbol=:zero)
    observable === :logical_x && return "Logical X error rate"
    observable === :logical_z && return "Logical Z error rate"
    if observable === :state_failure
        return "State-failure rate ($(_comparison_state_label(logical_state)))"
    end
    return "Logical failure rate ($(observable))"
end

function _plot_channel_curves!(axis, scan::ChannelFailureScan; scaled_fit=nothing)
    colors = (:steelblue, :darkorange, :seagreen, :purple, :firebrick)
    for (index, distance_value) in enumerate(scan.distances)
        points = filter(point -> point.distance == distance_value, scan.points)
        x = scaled_fit === nothing ?
            [_physical_error_rate(scan, point) for point in points] :
            [scaled_error_rate(_physical_error_rate(scan, point), point.distance, scaled_fit)
             for point in points]
        rates = channel_failure_rate.(Ref(scan), points)
        errors = channel_failure_standard_error.(Ref(scan), points)
        color = colors[mod1(index, length(colors))]
        CairoMakie.lines!(axis, x, rates; color, linewidth=2, label="d = $distance_value")
        CairoMakie.scatter!(axis, x, rates; color, markersize=8)
        CairoMakie.errorbars!(axis, x, rates, errors; color, whiskerwidth=7)
    end
end

function _comparison_panel!(position, scan::ChannelFailureScan, fit::ThresholdFit; title::AbstractString)
    suffix = fit.exploratory ? " — Exploratory" : ""
    axis = CairoMakie.Axis(
        position;
        title=String(title) * " — " * _comparison_state_label(scan.logical_state) * suffix,
        xlabel="Physical error rate",
        ylabel=_logical_axis_label(scan.logical_observable, scan.logical_state),
    )
    _plot_channel_curves!(axis, scan)
    CairoMakie.axislegend(axis; position=:lt)

    if fit.status === :success
        CairoMakie.vlines!(axis, [something(fit.p_c)]; color=:black, linestyle=:dash,
            linewidth=2)
        pc_error = fit.p_c_standard_error === nothing ? "" :
            " ± $(round(fit.p_c_standard_error; sigdigits=2))"
        nu_error = fit.nu_standard_error === nothing ? "" :
            " ± $(round(fit.nu_standard_error; sigdigits=2))"
        annotation = "p₍c₎ = $(round(something(fit.p_c); sigdigits=3))$(pc_error)\n" *
            "ν = $(round(something(fit.nu); sigdigits=3))$(nu_error)"
        CairoMakie.text!(axis, 0.97, 0.06; text=annotation, space=:relative,
            align=(:right, :bottom), fontsize=13)
        inset = CairoMakie.Axis(
            position;
            width=CairoMakie.Relative(0.43), height=CairoMakie.Relative(0.38),
            halign=0.97, valign=0.96,
            xlabel="(p − p₍c₎) d¹ᐟᵛ", ylabel=scan.logical_observable === :state_failure ?
                "state failure" : "logical error", titlesize=11,
            xlabelsize=9, ylabelsize=9, xticklabelsize=8, yticklabelsize=8,
            backgroundcolor=(:white, 0.92),
        )
        _plot_channel_curves!(inset, scan; scaled_fit=fit)
    else
        CairoMakie.text!(axis, 0.03, 0.05;
            text="Fit unavailable: $(fit.diagnostic)", space=:relative,
            align=(:left, :bottom), fontsize=12)
    end
    return axis
end

"""Render the 2×2 As/Bp × X-only/Z-only comparison figure."""
function plot_construction_channel_comparison(fitted::FittedConstructionChannelComparison)
    figure = CairoMakie.Figure(size=(1400, 980))
    constructions = (:as, :bp)
    channels = (:x_only, :z_only)
    for (row, construction) in enumerate(constructions), (column, error_channel) in enumerate(channels)
        scan = comparison_series(fitted.raw, construction, error_channel)
        fit = _comparison_fit(fitted, construction, error_channel)
        _comparison_panel!(
            figure[row, column], scan, fit;
            title="$(_construction_label(construction)) — $(_channel_label(error_channel))")
    end
    return figure
end

"""Render one validated construction/noise-channel comparison panel."""
function plot_construction_channel_panel(
        fitted::FittedConstructionChannelComparison,
        construction::Symbol,
        error_channel::Symbol)
    scan = comparison_series(fitted.raw, construction, error_channel)
    fit = _comparison_fit(fitted, construction, error_channel)
    figure = CairoMakie.Figure(size=(700, 560))
    _comparison_panel!(
        figure[1, 1], scan, fit;
        title="$(_construction_label(construction)) — $(_channel_label(error_channel))")
    return figure
end

function _comparison_basename(basename::AbstractString)
    isempty(basename) && throw(ArgumentError("basename must not be empty"))
    (basename in (".", "..") || occursin('/', basename) || occursin('\\', basename)) &&
        throw(ArgumentError("basename must be a filename without path components"))
    return String(basename)
end

function _save_comparison_figure(paths, figure)
    CairoMakie.save(paths.svg, figure)
    CairoMakie.save(paths.pdf, figure)
    CairoMakie.save(paths.png, figure)
    return paths
end

"""Save comparison CSVs, one combined figure, and all four standalone panels."""
function save_construction_channel_comparison(
        fitted::FittedConstructionChannelComparison,
        output_dir::AbstractString;
        basename::AbstractString="rotated_planar_comparison")
    leaf = _comparison_basename(basename)
    directory = abspath(output_dir)
    mkpath(directory)
    raw_csv = joinpath(directory, "$(leaf)-raw.csv")
    fits_csv = joinpath(directory, "$(leaf)-fits.csv")
    write_comparison_raw_csv(raw_csv, fitted)
    write_comparison_fit_csv(fits_csv, fitted)

    combined_paths = (
        svg=joinpath(directory, "$(leaf)-combined.svg"),
        pdf=joinpath(directory, "$(leaf)-combined.pdf"),
        png=joinpath(directory, "$(leaf)-combined.png"),
    )
    _save_comparison_figure(combined_paths, plot_construction_channel_comparison(fitted))

    panels = Dict{Tuple{Symbol,Symbol}, NamedTuple{(:svg, :pdf, :png),
        Tuple{String,String,String}}}()
    for (construction, error_channel) in _COMPARISON_SERIES_ORDER
        stem = "$(leaf)-$(construction)-$(replace(string(error_channel), '_' => '-'))"
        panel_paths = (
            svg=joinpath(directory, "$(stem).svg"),
            pdf=joinpath(directory, "$(stem).pdf"),
            png=joinpath(directory, "$(stem).png"),
        )
        panels[(construction, error_channel)] = _save_comparison_figure(
            panel_paths,
            plot_construction_channel_panel(fitted, construction, error_channel))
    end
    return (raw_csv=raw_csv, fits_csv=fits_csv, combined=combined_paths, panels=panels)
end
