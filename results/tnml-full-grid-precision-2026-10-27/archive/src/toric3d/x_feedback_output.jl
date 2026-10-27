const _TORIC3D_X_FEEDBACK_BASENAME = "toric3d-x-feedback"

function _toric3d_x_feedback_output_stem(directory::AbstractString,
        basename::AbstractString)
    isempty(strip(directory)) && throw(ArgumentError("output directory must not be empty"))
    !isempty(basename) && Base.basename(basename) == basename ||
        throw(ArgumentError("basename must be a nonempty filename stem"))
    absolute = abspath(directory)
    isfile(absolute) && throw(ArgumentError("output directory is an existing file"))
    mkpath(absolute)
    return joinpath(absolute, basename)
end

function _write_toric3d_x_feedback_raw(scan::Toric3DXFeedbackScan, path::AbstractString)
    columns = ["steps", "Ly", "Lz", "qubits", "p", "metric", "shots",
        "seed", "point_seed", "routing", "final_only_failures",
        "slice_feedback_failures", "both_failures", "final_only_only",
        "slice_feedback_only", "final_only_rate", "final_only_ci_low",
        "final_only_ci_high", "slice_feedback_rate", "slice_feedback_ci_low",
        "slice_feedback_ci_high", "difference_feedback_minus_final",
        "difference_standard_error", "difference_ci_low", "difference_ci_high",
        "elapsed_seconds", "ldpc_version", "max_iter", "bp_method",
        "ms_scaling_factor", "schedule", "omp_thread_count", "osd_method",
        "osd_order"]
    settings = scan.settings
    _toric3d_atomic_write(path) do io
        println(io, join(columns, ','))
        for point in scan.points
            steps, Ly, Lz = point.dimensions
            for (metric_name, metric) in pairs(point.metrics)
                values = (steps, Ly, Lz, 3 * (steps + 1) * Ly * Lz, point.p,
                    metric_name, point.shots, something(scan.seed, "unspecified"),
                    _toric3d_x_feedback_seed_string(point.point_seed), point.routing,
                    metric.final_only_failures, metric.layer_feedback_failures,
                    metric.both_failures, metric.final_only_only,
                    metric.layer_feedback_only, metric.final_only_rate,
                    metric.final_only_ci..., metric.layer_feedback_rate,
                    metric.layer_feedback_ci..., metric.difference,
                    metric.difference_standard_error, metric.difference_ci...,
                    point.elapsed_seconds, point.decoder_version,
                    settings.max_iter, settings.bp_method,
                    settings.ms_scaling_factor, settings.schedule,
                    settings.omp_thread_count, settings.osd_method,
                    settings.osd_order)
                println(io, join(_toric3d_csv.(values), ','))
            end
        end
    end
    return path
end

function _toric3d_x_feedback_complete(scan::Toric3DXFeedbackScan)
    expected = Set((dimensions, p) for dimensions in scan.dimensions
        for p in scan.error_rates)
    observed = Set((point.dimensions, point.p) for point in scan.points)
    return length(scan.points) == length(expected) && observed == expected
end

function _toric3d_x_feedback_metadata(scan::Toric3DXFeedbackScan)
    expected_points = length(scan.dimensions) * length(scan.error_rates)
    versions = sort!(unique(point.decoder_version for point in scan.points))
    return Dict{String,Any}(
        "format_version" => 1,
        "dimensions_steps_Ly_Lz" => [collect(dimensions) for dimensions in scan.dimensions],
        "error_rates" => scan.error_rates,
        "shots_per_point" => scan.shots,
        "seed_metadata" => something(scan.seed, "unspecified"),
        "rng_type" => scan.rng_type,
        "point_rng_type" => "Random.Xoshiro",
        "point_seeds_hex" => _toric3d_x_feedback_seed_string.(scan.point_seeds),
        "started_at_unix" => scan.started_at_unix,
        "completed_points" => length(scan.points),
        "expected_points" => expected_points,
        "complete" => _toric3d_x_feedback_complete(scan),
        "batch_size" => scan.batch_size,
        "routing" => string(scan.routing),
        "resumed" => scan.resumed,
        "resume_supported" => true,
        "checkpoint_filename" => _TORIC3D_X_FEEDBACK_CHECKPOINT,
        "protocols" => ["final_only", "slice_feedback"],
        "logical_qubits" => 2,
        "initial_state" => "logical |00> prepared by ideal X-star measurement and ideal Z correction",
        "state_failure" => "either residual logical-X membrane parity after final recovery",
        "logical_diagnostics" => ["logical_1", "logical_2", "state_failure"],
        "geometry" => "diagonal encoded slice with n=3(steps+1)LyLz and two logical qubits",
        "transfer_schedule" => "seven-gate diagonal transfer: H followed by six CNOTs",
        "active_gate_outputs_per_transfer_cell" => 13,
        "seed_gate_fault_locations" => 0,
        "noise" => "independent Bernoulli X fault of rate p after each participating H/CNOT output",
        "ideal_operations" => ["initialization", "idle qubits", "syndrome measurements",
            "classical processing", "Pauli corrections"],
        "feedback" => "ideal active-slice star measurement and MWPM correction after each diagonal slice",
        "final_decoder" => "full-plaquette BP+OSD X recovery, validated against the measured syndrome",
        "decoder_prior" => "uniform declared physical p on every final data qubit",
        "decoder_limitations" => "uniform priors neglect nonuniform exposure across qubits and circuit/feedback correlations",
        "ldpc_versions" => versions,
        "individual_interval" => "95% Wilson score interval",
        "paired_difference" => "slice_feedback minus final_only using paired discordances",
        "paired_interval" => "conservative approximate 95% interval from 97.5% Wilson discordance intervals",
        "decoder_settings" => Dict(string(key) => value
            for (key, value) in pairs(_toric3d_settings_metadata(scan.settings))))
end

_toric3d_x_feedback_size_label(dimensions) = all(==(dimensions[1]), dimensions) ?
    "L = $(dimensions[1])" : join(dimensions, " × ")

function _toric3d_x_feedback_series(scan, dimensions)
    return sort(filter(point -> point.dimensions == dimensions, scan.points);
        by=point -> point.p)
end

"""Plot state failure, individual logical rates, and paired protocol differences."""
function plot_toric3d_x_feedback_threshold(scan::Toric3DXFeedbackScan)::CairoMakie.Figure
    isempty(scan.points) && throw(ArgumentError("scan has no completed points"))
    figure = CairoMakie.Figure(size=(1500, 570))
    colors = ["#2364AA", "#D1495B", "#258477", "#8B5FBF", "#A66326",
        "#5F6B6D"]
    state_axis = CairoMakie.Axis(figure[1, 1];
        xlabel="X-fault probability p per gate output",
        ylabel="Logical |00> state-failure probability",
        title="Either logical-X membrane fails")
    logical_axis = CairoMakie.Axis(figure[1, 2];
        xlabel="X-fault probability p per gate output",
        ylabel="Residual logical-X probability",
        title="Individual logicals (● logical 1, ◆ logical 2)")
    difference_axis = CairoMakie.Axis(figure[1, 3];
        xlabel="X-fault probability p per gate output",
        ylabel="Feedback − final-only",
        title="Paired state-failure difference")
    CairoMakie.hlines!(difference_axis, [0.0]; color=(:gray40, 0.8),
        linestyle=:dot, linewidth=1.5)

    size_handles = Any[]
    size_labels = String[]
    for (size_index, dimensions) in enumerate(scan.dimensions)
        points = _toric3d_x_feedback_series(scan, dimensions)
        isempty(points) && continue
        color = colors[mod1(size_index, length(colors))]
        x = [point.p for point in points]
        state = [point.metrics.state_failure for point in points]
        for (protocol, linestyle, marker) in
                ((:final_only, :solid, :circle), (:layer_feedback, :dash, :rect))
            y = [getproperty(metric, Symbol(protocol, :_rate)) for metric in state]
            intervals = [getproperty(metric, Symbol(protocol, :_ci)) for metric in state]
            CairoMakie.errorbars!(state_axis, x, y,
                [y[index] - intervals[index][1] for index in eachindex(y)],
                [intervals[index][2] - y[index] for index in eachindex(y)];
                color=(color, 0.55), whiskerwidth=4)
            CairoMakie.lines!(state_axis, x, y; color, linestyle, linewidth=2)
            CairoMakie.scatter!(state_axis, x, y; color, marker, markersize=7)
        end

        for (metric_name, marker) in ((:logical_1, :circle), (:logical_2, :diamond))
            metrics = [getproperty(point.metrics, metric_name) for point in points]
            for (protocol, linestyle) in
                    ((:final_only, :solid), (:layer_feedback, :dash))
                y = [getproperty(metric, Symbol(protocol, :_rate)) for metric in metrics]
                intervals = [getproperty(metric, Symbol(protocol, :_ci)) for metric in metrics]
                CairoMakie.errorbars!(logical_axis, x, y,
                    [y[index] - intervals[index][1] for index in eachindex(y)],
                    [intervals[index][2] - y[index] for index in eachindex(y)];
                    color=(color, 0.4), whiskerwidth=3)
                CairoMakie.lines!(logical_axis, x, y; color, linestyle, linewidth=1.8)
                CairoMakie.scatter!(logical_axis, x, y; color, marker, markersize=6)
            end
        end

        difference = [metric.difference for metric in state]
        difference_ci = [metric.difference_ci for metric in state]
        CairoMakie.errorbars!(difference_axis, x, difference,
            [difference[index] - difference_ci[index][1] for index in eachindex(difference)],
            [difference_ci[index][2] - difference[index] for index in eachindex(difference)];
            color=(color, 0.65), whiskerwidth=4)
        CairoMakie.lines!(difference_axis, x, difference; color, linewidth=2)
        CairoMakie.scatter!(difference_axis, x, difference; color,
            marker=:circle, markersize=7)
        push!(size_handles, CairoMakie.LineElement(color=color, linewidth=3))
        push!(size_labels, _toric3d_x_feedback_size_label(dimensions))
    end

    handles = vcat(size_handles, Any[
        CairoMakie.LineElement(color=:gray25, linestyle=:solid, linewidth=2),
        CairoMakie.LineElement(color=:gray25, linestyle=:dash, linewidth=2)])
    labels = vcat(size_labels, ["Final-only", "Slice feedback"])
    CairoMakie.Legend(figure[2, 1:3], handles, labels;
        nbanks=min(length(labels), 4), framevisible=false, tellheight=true)
    CairoMakie.Label(figure[3, 1:3],
        "$(scan.shots) paired shots/point · routing=$(scan.routing) · ideal measurement and correction · 95% intervals\n" *
        "Final BP+OSD uses uniform p priors and neglects nonuniform exposure and circuit/feedback correlations.";
        fontsize=12, tellwidth=false)
    return figure
end

function _toric3d_x_feedback_save_figure(figure, path::AbstractString)
    extension = splitext(path)[2]
    temporary = path * ".tmp" * extension
    CairoMakie.save(temporary, figure)
    mv(temporary, path; force=true)
    return path
end

"""Save the raw CSV, metadata, and PDF/SVG/PNG diagonal-feedback figures."""
function save_toric3d_x_feedback_threshold(scan::Toric3DXFeedbackScan,
        directory::AbstractString;
        basename::AbstractString=_TORIC3D_X_FEEDBACK_BASENAME)::NamedTuple
    stem = _toric3d_x_feedback_output_stem(directory, basename)
    csv = _write_toric3d_x_feedback_raw(scan, stem * "-raw.csv")
    metadata = stem * "-metadata.toml"
    _toric3d_atomic_write(metadata) do io
        TOML.print(io, _toric3d_x_feedback_metadata(scan); sorted=true)
    end
    figure = plot_toric3d_x_feedback_threshold(scan)
    svg = _toric3d_x_feedback_save_figure(figure, stem * ".svg")
    pdf = _toric3d_x_feedback_save_figure(figure, stem * ".pdf")
    png = _toric3d_x_feedback_save_figure(figure, stem * ".png")
    return (csv=csv, metadata=metadata, svg=svg, pdf=pdf, png=png)
end
