"""A deterministic rectangular scan of circuit-level logical failure rates."""
struct LogicalFailureScan
    distances::Vector{Int}
    error_rates::Vector{Float64}
    logical_state::Symbol
    construction::Symbol
    boundary_orientation::Symbol
    clock::Symbol
    shots::Int
    batch_size::Int
    seed::Union{Nothing,Int}
    points::Vector{LogicalFailurePoint}
end

function _scan_distances(values)
    distances = try
        collect(values)
    catch
        throw(ArgumentError("distances must be a nonempty iterable of integers"))
    end
    isempty(distances) && throw(ArgumentError("distances must not be empty"))
    all(value -> value isa Integer, distances) ||
        throw(ArgumentError("distances must contain only integers"))
    converted = try
        Int.(distances)
    catch
        throw(ArgumentError("distances must be representable as Int"))
    end
    all(value -> value >= 3 && isodd(value), converted) ||
        throw(ArgumentError("distances must be odd integers at least 3"))
    all(diff(converted) .> 0) ||
        throw(ArgumentError("distances must be unique and strictly increasing"))
    return converted
end

function _scan_error_rates(values)
    error_rates = try
        collect(values)
    catch
        throw(ArgumentError("error_rates must be a nonempty iterable of real numbers"))
    end
    isempty(error_rates) && throw(ArgumentError("error_rates must not be empty"))
    all(value -> value isa Real, error_rates) ||
        throw(ArgumentError("error_rates must contain only real numbers"))
    converted = Float64.(error_rates)
    all(value -> isfinite(value) && 0 <= value < 0.5, converted) ||
        throw(ArgumentError(
            "error_rates must be finite and satisfy 0 <= error rate < 0.5"))
    all(diff(converted) .> 0) ||
        throw(ArgumentError("error_rates must be unique and strictly increasing"))
    return converted
end

"""
    scan_logical_failure(rng; distances, error_rates, kwargs...) -> LogicalFailureScan

Estimate every `(distance, error_rate)` pair in distance-major order. The
supplied RNG is consumed continuously across the scan; `seed` is metadata only.
By default, As prepares `|+_L>` and Bp prepares `|0_L>` (`logical_state=:native`).
Their ideal product inputs are all-plus and all-zero respectively.
"""
function scan_logical_failure(
        rng::Random.AbstractRNG;
        distances,
        error_rates,
        logical_state::Symbol=:native,
        construction::Symbol=:bp,
        boundary_orientation::Symbol=:x_ns,
        clock::Symbol=:gate_layer,
        shots::Integer=10_000,
        batch_size::Integer=10_000,
        seed=nothing,
        progress_io=nothing,
    )::LogicalFailureScan
    distance_values = _scan_distances(distances)
    error_rate_values = _scan_error_rates(error_rates)
    shot_count = _positive_machine_int(shots, "shots")
    chunk_limit = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)
    (progress_io === nothing || progress_io isa IO) ||
        throw(ArgumentError("progress_io must be nothing or an IO stream"))

    validation_code = RotatedPlanarCode(
        first(distance_values); boundary_orientation)
    logical_state = rotated_planar_encoder(
        validation_code; logical_state, construction).logical_state
    CircuitPauliNoise(first(error_rate_values); clock)

    points = LogicalFailurePoint[]
    sizehint!(points, length(distance_values) * length(error_rate_values))
    for distance_value in distance_values
        code = RotatedPlanarCode(distance_value; boundary_orientation)
        encoder = rotated_planar_encoder(code; logical_state, construction)
        for error_rate in error_rate_values
            noise = CircuitPauliNoise(error_rate; clock)
            point = estimate_logical_failure(
                rng, encoder, noise;
                shots=shot_count, batch_size=chunk_limit, seed=seed_value)
            push!(points, point)
            progress_io === nothing || println(
                progress_io,
                "completed d=$(point.distance) p=$(error_rate) " *
                "any=$(point.any_logical_failure_rate)")
        end
    end

    return LogicalFailureScan(
        distance_values, error_rate_values, logical_state, construction,
        boundary_orientation, clock, shot_count, chunk_limit, seed_value, points)
end

const _LOGICAL_FAILURE_CSV_COLUMNS = (
    :distance,
    :logical_state,
    :construction,
    :boundary_orientation,
    :clock,
    :p_x,
    :p_z,
    :shots,
    :seed,
    :logical_x_failures,
    :logical_x_failure_rate,
    :logical_x_standard_error,
    :logical_z_failures,
    :logical_z_failure_rate,
    :logical_z_standard_error,
    :any_logical_failures,
    :any_logical_failure_rate,
    :any_logical_standard_error,
    :state_failures,
    :state_failure_rate,
    :state_standard_error,
)

_csv_value(value) = value === nothing ? "" : string(value)

"""Write a logical-failure scan as dependency-free tidy CSV."""
function write_logical_failure_csv(
        path::AbstractString, scan::LogicalFailureScan)
    output_path = abspath(path)
    open(output_path, "w") do io
        println(io, join(string.(_LOGICAL_FAILURE_CSV_COLUMNS), ','))
        for point in scan.points
            println(io, join(
                (_csv_value(getfield(point, column))
                 for column in _LOGICAL_FAILURE_CSV_COLUMNS),
                ','))
        end
    end
    return output_path
end

"""Plot logical-X, logical-Z, and either-logical scan results with errors."""
function plot_logical_failure_scan(scan::LogicalFailureScan)
    figure = CairoMakie.Figure(size=(1320, 420))
    panels = (
        (:logical_x_failure_rate, :logical_x_standard_error,
         "Logical-X failure rate"),
        (:logical_z_failure_rate, :logical_z_standard_error,
         "Logical-Z failure rate"),
        (:any_logical_failure_rate, :any_logical_standard_error,
         "Either-logical failure rate"),
    )
    colors = (:steelblue, :darkorange, :seagreen, :purple, :firebrick)

    for (panel_index, (rate_field, error_field, title)) in enumerate(panels)
        axis = CairoMakie.Axis(
            figure[1, panel_index];
            title,
            xlabel="Physical error rate",
            ylabel="Failure rate",
        )
        for (distance_index, distance_value) in enumerate(scan.distances)
            points = filter(point -> point.distance == distance_value, scan.points)
            x = [point.p_x for point in points]
            rates = [getfield(point, rate_field) for point in points]
            errors = [getfield(point, error_field) for point in points]
            color = colors[mod1(distance_index, length(colors))]
            CairoMakie.lines!(
                axis, x, rates; color, linewidth=2, label="d = $distance_value")
            CairoMakie.scatter!(axis, x, rates; color, markersize=9)
            CairoMakie.errorbars!(axis, x, rates, errors; color, whiskerwidth=8)
        end
        CairoMakie.axislegend(axis; position=:lt)
    end
    return figure
end

"""Save CSV, SVG, PDF, and PNG artifacts for a logical-failure scan."""
function save_logical_failure_scan(
        scan::LogicalFailureScan,
        output_dir::AbstractString;
        basename::AbstractString="rotated_planar_logical_failure")
    isempty(basename) && throw(ArgumentError("basename must not be empty"))
    (basename in (".", "..") || occursin('/', basename) || occursin('\\', basename)) &&
        throw(ArgumentError("basename must be a filename without path components"))
    directory = abspath(output_dir)
    mkpath(directory)
    paths = (
        csv=joinpath(directory, "$basename.csv"),
        svg=joinpath(directory, "$basename.svg"),
        pdf=joinpath(directory, "$basename.pdf"),
        png=joinpath(directory, "$basename.png"),
    )
    write_logical_failure_csv(paths.csv, scan)
    figure = plot_logical_failure_scan(scan)
    CairoMakie.save(paths.svg, figure)
    CairoMakie.save(paths.pdf, figure)
    CairoMakie.save(paths.png, figure)
    return paths
end
