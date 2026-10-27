const _TORIC3D_DIAGONAL_FEEDBACK_BASENAME = "toric3d-diagonal-feedback"

const _TORIC3D_DIAGONAL_FEEDBACK_COLUMNS = [
    "steps", "Ly", "Lz", "qubits", "d_X", "d_Z", "p", "channel",
    "metric", "shots", "seed", "point_seed", "x_routing",
    "joint_n0_s0_f0", "joint_n1_s0_f0", "joint_n0_s1_f0",
    "joint_n1_s1_f0", "joint_n0_s0_f1", "joint_n1_s0_f1",
    "joint_n0_s1_f1", "joint_n1_s1_f1",
    "none_failures", "none_rate", "none_ci_low", "none_ci_high",
    "slice_local_failures", "slice_local_rate", "slice_local_ci_low",
    "slice_local_ci_high", "full_prefix_failures", "full_prefix_rate",
    "full_prefix_ci_low", "full_prefix_ci_high",
    "none_vs_slice_local_neither", "none_vs_slice_local_none_only",
    "none_vs_slice_local_slice_local_only", "none_vs_slice_local_both",
    "none_vs_slice_local_difference_slice_local_minus_none",
    "none_vs_slice_local_standard_error", "none_vs_slice_local_ci_low",
    "none_vs_slice_local_ci_high", "none_vs_full_prefix_neither",
    "none_vs_full_prefix_none_only", "none_vs_full_prefix_full_prefix_only",
    "none_vs_full_prefix_both",
    "none_vs_full_prefix_difference_full_prefix_minus_none",
    "none_vs_full_prefix_standard_error", "none_vs_full_prefix_ci_low",
    "none_vs_full_prefix_ci_high", "slice_local_vs_full_prefix_neither",
    "slice_local_vs_full_prefix_slice_local_only",
    "slice_local_vs_full_prefix_full_prefix_only",
    "slice_local_vs_full_prefix_both",
    "slice_local_vs_full_prefix_difference_full_prefix_minus_slice_local",
    "slice_local_vs_full_prefix_standard_error",
    "slice_local_vs_full_prefix_ci_low", "slice_local_vs_full_prefix_ci_high",
    "decoder_build_seconds", "sampling_seconds", "replay_seconds",
    "pymatching_version", "ldpc_version", "invalid_recoveries",
    "max_iter", "bp_method", "ms_scaling_factor", "schedule",
    "omp_thread_count", "osd_method", "osd_order"]

function _toric3d_diagonal_feedback_output_stem(directory::AbstractString,
        basename::AbstractString)
    isempty(strip(directory)) &&
        throw(ArgumentError("output directory must not be empty"))
    !isempty(basename) && Base.basename(basename) == basename ||
        throw(ArgumentError("basename must be a nonempty filename stem"))
    absolute = abspath(directory)
    isfile(absolute) && throw(ArgumentError("output directory is an existing file"))
    mkpath(absolute)
    return joinpath(absolute, basename)
end

_toric3d_protocol_values(estimate::Toric3DProtocolEstimate) =
    (estimate.failures, estimate.rate, estimate.ci...)

_toric3d_pair_values(estimate::Toric3DPairEstimate) =
    (estimate.neither, estimate.first_only, estimate.second_only, estimate.both,
        estimate.difference, estimate.standard_error, estimate.ci...)

function _write_toric3d_diagonal_feedback_raw(
        scan::Toric3DDiagonalFeedbackScan, path::AbstractString)
    settings = scan.settings
    _toric3d_atomic_write(path) do io
        println(io, join(_TORIC3D_DIAGONAL_FEEDBACK_COLUMNS, ','))
        for point in scan.points
            steps, Ly, Lz = point.dimensions
            for (metric_name, metric) in pairs(point.metrics)
                values = (steps, Ly, Lz, 3 * (steps + 1) * Ly * Lz,
                    (steps + 1) * min(Ly, Lz), 2 * min(Ly, Lz), point.p,
                    point.channel, metric_name, point.shots,
                    something(scan.seed, "unspecified"),
                    _toric3d_feedback_seed_string(point.point_seed), scan.x_routing,
                    metric.joint_counts...,
                    _toric3d_protocol_values(metric.none)...,
                    _toric3d_protocol_values(metric.slice_local)...,
                    _toric3d_protocol_values(metric.full_prefix)...,
                    _toric3d_pair_values(metric.none_vs_slice_local)...,
                    _toric3d_pair_values(metric.none_vs_full_prefix)...,
                    _toric3d_pair_values(metric.slice_local_vs_full_prefix)...,
                    point.phase_seconds.decoder_build,
                    point.phase_seconds.sampling, point.phase_seconds.replay,
                    point.decoder_versions.pymatching,
                    point.decoder_versions.ldpc, point.invalid_recoveries,
                    settings.max_iter, settings.bp_method,
                    settings.ms_scaling_factor, settings.schedule,
                    settings.omp_thread_count, settings.osd_method,
                    settings.osd_order)
                length(values) == length(_TORIC3D_DIAGONAL_FEEDBACK_COLUMNS) ||
                    error("diagonal-feedback CSV schema/value mismatch")
                println(io, join(_toric3d_csv.(values), ','))
            end
        end
    end
    return path
end

function _toric3d_diagonal_feedback_complete(scan::Toric3DDiagonalFeedbackScan)
    expected = Set((dimensions, p, channel) for dimensions in scan.dimensions
        for p in scan.error_rates for channel in scan.channels)
    observed = Set((point.dimensions, point.p, point.channel)
        for point in scan.points)
    return length(scan.points) == length(expected) && observed == expected
end

function _toric3d_fnv1a64(path::AbstractString)
    value = UInt64(0xcbf29ce484222325)
    for byte in read(path)
        value = (value ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    end
    return string(value; base=16, pad=16)
end

function _toric3d_diagonal_feedback_provenance()
    root = normpath(joinpath(@__DIR__, "..", ".."))
    sources = [
        "src/toric3d/feedback_strength_geometry.jl",
        "src/toric3d/feedback_strength_decoders.jl",
        "src/toric3d/feedback_strength_x.jl",
        "src/toric3d/feedback_strength_z.jl",
        "src/toric3d/feedback_strength_z_batch.jl",
        "src/toric3d/feedback_strength_scan.jl",
        "src/toric3d/feedback_strength_output.jl",
        "examples/compare_toric3d_diagonal_feedback.jl"]
    hashes = Dict{String,String}()
    for relative in sources
        path = joinpath(root, relative)
        isfile(path) && (hashes[relative] = _toric3d_fnv1a64(path))
    end
    git(args...) = readchomp(Cmd(vcat(["git", "-C", root], collect(args))))
    state = try
        status = git("status", "--porcelain")
        Dict{String,Any}(
            "commit" => git("rev-parse", "HEAD"),
            "branch" => git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty" => !isempty(status))
    catch
        Dict{String,Any}(
            "commit" => "unavailable", "branch" => "unavailable",
            "dirty" => false)
    end
    return hashes, state
end

function _toric3d_diagonal_feedback_fault_locations(scan)
    entries = Dict{String,Any}[]
    for dimensions in scan.dimensions, channel in scan.channels
        steps, Ly, Lz = dimensions
        transfer = 13 * steps * Ly * Lz
        seed = channel === :z_only ? 3 * Ly * Lz : 0
        push!(entries, Dict{String,Any}(
            "dimensions_steps_Ly_Lz" => collect(dimensions),
            "channel" => string(channel), "seed_gate_outputs" => seed,
            "transfer_gate_outputs" => transfer,
            "total_gate_outputs" => seed + transfer))
    end
    return entries
end

function _toric3d_diagonal_feedback_channel_models(scan)
    return Dict{String,Any}(
        "x_only" => Dict{String,Any}(
            "initial_state" => "logical |00> from ideal seed X-star measurement and ideal Z correction",
            "fault_pauli" => "X",
            "seed_circuit" => "ideal encoded seed; no noisy seed gates",
            "seed_gate_outputs" => 0,
            "seed_gate_outputs_formula" => "0",
            "state_failure" => "either residual logical-X membrane",
            "slice_local_check_set" => "active-slice Z compatibility checks before each transfer",
            "slice_local_correction_support" =>
                "fresh transfer targets reached by the selected target feed-forward route",
            "slice_local_policy" =>
                "one-step target feed-forward; x_routing=full_path is a sensitivity only",
            "full_prefix_policy" =>
                "all completed prefix plaquettes, recovery on all mapped prefix qubits",
            "final_decoder" => "validated full-code BP+OSD X recovery",
            "decoder_prior" => "uniform declared p, clamped to [1e-12,1-1e-12]",
            "decoder_limitations" =>
                "neglects nonuniform exposure, circuit correlations, and feedback-conditioned correlations"),
        "z_only" => Dict{String,Any}(
            "initial_state" => "logical |++> from noisy H on every seed edge and ideal compatibility correction",
            "fault_pauli" => "Z",
            "seed_circuit" => "one H output on every one of 3*Ly*Lz seed edges",
            "seed_gate_outputs" => "dimension dependent",
            "seed_gate_outputs_formula" => "3*Ly*Lz",
            "state_failure" => "either residual logical-Z winding",
            "slice_local_check_set" =>
                "seed checks, then a three-check-sheet moving window",
            "slice_local_correction_support" =>
                "seed edges at r=0; newest two data slices r-1:r thereafter",
            "slice_local_policy" => "complete deterministic unit-weight MWPM window recovery",
            "full_prefix_policy" =>
                "all completed prefix stars, recovery on all mapped prefix qubits",
            "final_decoder" => "validated full-code deterministic unit-weight MWPM Z recovery",
            "decoder_prior" => "unit data-edge weights",
            "decoder_limitations" =>
                "neglects circuit correlations and feedback-conditioned correlations"))
end

function _toric3d_diagonal_feedback_metadata(scan::Toric3DDiagonalFeedbackScan)
    expected = length(scan.dimensions) * length(scan.error_rates) *
        length(scan.channels)
    hashes, git_state = _toric3d_diagonal_feedback_provenance()
    pymatching = sort!(unique(point.decoder_versions.pymatching
        for point in scan.points))
    ldpc = sort!(unique(point.decoder_versions.ldpc for point in scan.points))
    return Dict{String,Any}(
        "format_version" => 1,
        "dimensions_steps_Ly_Lz" => [collect(value) for value in scan.dimensions],
        "error_rates" => scan.error_rates, "channels" => string.(scan.channels),
        "shots_per_point" => scan.shots, "batch_size" => scan.batch_size,
        "seed_metadata" => something(scan.seed, "unspecified"),
        "rng_type" => scan.rng_type, "point_rng_type" => "Random.Xoshiro",
        "rng_policy" =>
            "one deterministic point substream; each immutable fault record is replayed through all protocols",
        "point_seeds_hex" => _toric3d_feedback_seed_string.(scan.point_seeds),
        "started_at_unix" => scan.started_at_unix,
        "x_routing" => string(scan.x_routing), "resumed" => scan.resumed,
        "resume_supported" => true,
        "checkpoint_filename" => _TORIC3D_DIAGONAL_FEEDBACK_CHECKPOINT,
        "failure_filename" => _TORIC3D_DIAGONAL_FEEDBACK_FAILURE,
        "completed_points" => length(scan.points), "expected_points" => expected,
        "complete" => _toric3d_diagonal_feedback_complete(scan),
        "invalid_recoveries" => scan.invalid_recoveries,
        "decoder_builds_x" => scan.decoder_builds.x,
        "decoder_builds_z" => scan.decoder_builds.z,
        "protocols" => ["none", "slice_local", "full_prefix"],
        "metrics" => string.(collect(_TORIC3D_FEEDBACK_METRICS)),
        "logical_qubits" => 2,
        "geometry" =>
            "diagonal encoded slice with n=3(steps+1)LyLz, k=2, dX=(steps+1)min(Ly,Lz), dZ=2min(Ly,Lz)",
        "transfer_schedule" =>
            "seven gates per cell: H on first fresh edge, two fresh-edge CNOTs, four old-control CNOTs",
        "active_gate_outputs_per_transfer_cell" => 13,
        "ideal_operations" => ["initialization except declared Z seed H gates",
            "idle qubits", "syndrome measurements", "classical processing",
            "Pauli corrections"],
        "paired_replay" =>
            "the identical sampled preparation fault record is replayed through all three protocols",
        "joint_pattern_encoding" =>
            "mask = none + 2*slice_local + 4*full_prefix",
        "individual_interval" => "95% Wilson score interval",
        "paired_difference_convention" => "second protocol minus first protocol",
        "paired_interval" =>
            "conservative approximate 95% interval from 97.5% Wilson discordance intervals",
        "timing_fields" =>
            "aggregate decoder construction, fault sampling, and three-protocol replay wall seconds per point",
        "resource_scaling" => Dict(
            "slice_local" => "O(steps*Ly*Lz) check participation",
            "full_prefix" => "O(steps^2*Ly*Lz) check participation"),
        "channel_comparison_warning" =>
            "X and Z use different logical input states, seed circuits, fault locations, and decoders; compare protocols within a channel",
        "decoder_settings" => Dict(string(key) => value
            for (key, value) in pairs(_toric3d_settings_metadata(scan.settings))),
        "decoder_versions" => Dict("pymatching" => pymatching, "ldpc" => ldpc),
        "channel_models" => _toric3d_diagonal_feedback_channel_models(scan),
        "fault_locations" => _toric3d_diagonal_feedback_fault_locations(scan),
        "julia_version" => string(VERSION),
        "source_hash_algorithm" => "FNV-1a-64",
        "source_hashes" => hashes, "git_state" => git_state)
end

_toric3d_diagonal_size_label(dimensions) = all(==(dimensions[1]), dimensions) ?
    "L = $(dimensions[1])" : join(dimensions, " × ")

function _toric3d_diagonal_feedback_series(scan, dimensions, channel)
    return sort(filter(point -> point.dimensions == dimensions &&
        point.channel === channel, scan.points); by=point -> point.p)
end

function _toric3d_diagonal_protocol_style(protocol)
    protocol === :none && return (:solid, :circle,
        "No feedback (final decode only)")
    protocol === :slice_local && return (:dash, :rect, "Slice-local feedback")
    return (:dot, :utriangle, "Full-prefix feedback")
end

function _toric3d_plot_interval!(axis, x, estimates, color, linestyle, marker;
        label="")
    y = [estimate.rate for estimate in estimates]
    ci = [estimate.ci for estimate in estimates]
    CairoMakie.errorbars!(axis, x, y,
        [y[index] - ci[index][1] for index in eachindex(y)],
        [ci[index][2] - y[index] for index in eachindex(y)];
        color=(color, 0.38), whiskerwidth=3)
    CairoMakie.lines!(axis, x, y; color, linestyle, linewidth=2, label)
    CairoMakie.scatter!(axis, x, y; color, marker, markersize=7)
    return nothing
end

"""Plot three feedback strengths for each selected pure-Pauli noise channel."""
function plot_toric3d_diagonal_feedback(
        scan::Toric3DDiagonalFeedbackScan)::CairoMakie.Figure
    isempty(scan.points) && throw(ArgumentError("scan has no completed points"))
    channels = [channel for channel in scan.channels if
        any(point -> point.channel === channel, scan.points)]
    isempty(channels) && throw(ArgumentError("scan has no selected channel points"))
    colors = ["#2364AA", "#D1495B", "#258477", "#8B5FBF", "#A66326",
        "#5F6B6D"]
    figure = CairoMakie.Figure(size=(1550, 500 * length(channels) + 190))
    pair_styles = ((:none_vs_slice_local, :solid, :circle,
            "slice-local − none"),
        (:none_vs_full_prefix, :dash, :rect, "full-prefix − none"),
        (:slice_local_vs_full_prefix, :dot, :utriangle,
            "full-prefix − slice-local"))
    size_handles = Any[]
    size_labels = String[]

    for (row, channel) in enumerate(channels)
        pauli = channel === :x_only ? "X" : "Z"
        state = channel === :x_only ? "|00⟩" : "|++⟩"
        state_axis = CairoMakie.Axis(figure[row, 1];
            xlabel="$pauli-fault probability p per participating gate output",
            ylabel="Logical $state state-failure probability",
            title="$pauli-only · State failure")
        logical_axis = CairoMakie.Axis(figure[row, 2];
            xlabel="$pauli-fault probability p per participating gate output",
            ylabel="Residual logical-$pauli probability",
            title="$pauli-only · Individual logicals (● logical 1, ◆ logical 2)")
        difference_axis = CairoMakie.Axis(figure[row, 3];
            xlabel="$pauli-fault probability p per participating gate output",
            ylabel="Paired rate difference",
            title="$pauli-only · Paired differences (second − first)")
        CairoMakie.hlines!(difference_axis, [0.0]; color=(:gray40, 0.8),
            linestyle=:dashdot, linewidth=1.2)

        for (size_index, dimensions) in enumerate(scan.dimensions)
            points = _toric3d_diagonal_feedback_series(scan, dimensions, channel)
            isempty(points) && continue
            color = colors[mod1(size_index, length(colors))]
            x = [point.p for point in points]
            state_metrics = [point.metrics.state_failure for point in points]
            for protocol in (:none, :slice_local, :full_prefix)
                linestyle, marker, protocol_label =
                    _toric3d_diagonal_protocol_style(protocol)
                estimates = [getproperty(metric, protocol) for metric in state_metrics]
                _toric3d_plot_interval!(state_axis, x, estimates, color,
                    linestyle, marker; label=protocol_label)
            end

            for (metric_name, logical_marker) in
                    ((:logical_1, :circle), (:logical_2, :diamond))
                metrics = [getproperty(point.metrics, metric_name) for point in points]
                for protocol in (:none, :slice_local, :full_prefix)
                    linestyle, _, protocol_label =
                        _toric3d_diagonal_protocol_style(protocol)
                    estimates = [getproperty(metric, protocol) for metric in metrics]
                    _toric3d_plot_interval!(logical_axis, x, estimates, color,
                        linestyle, logical_marker;
                        label="$(protocol_label) · $(replace(string(metric_name), '_' => ' '))")
                end
            end

            for (pair_name, linestyle, marker, pair_label) in pair_styles
                pairs = [getproperty(metric, pair_name) for metric in state_metrics]
                y = [pair.difference for pair in pairs]
                ci = [pair.ci for pair in pairs]
                CairoMakie.errorbars!(difference_axis, x, y,
                    [y[index] - ci[index][1] for index in eachindex(y)],
                    [ci[index][2] - y[index] for index in eachindex(y)];
                    color=(color, 0.38), whiskerwidth=3)
                CairoMakie.lines!(difference_axis, x, y; color, linestyle,
                    linewidth=2, label=pair_label)
                CairoMakie.scatter!(difference_axis, x, y; color, marker,
                    markersize=7)
            end
            row == 1 && begin
                push!(size_handles, CairoMakie.LineElement(color=color, linewidth=3))
                push!(size_labels, _toric3d_diagonal_size_label(dimensions))
            end
        end
    end

    protocol_handles = Any[
        CairoMakie.LineElement(color=:gray25, linestyle=:solid, linewidth=2),
        CairoMakie.LineElement(color=:gray25, linestyle=:dash, linewidth=2),
        CairoMakie.LineElement(color=:gray25, linestyle=:dot, linewidth=2)]
    protocol_labels = ["No feedback (final decode only)",
        "Slice-local feedback", "Full-prefix feedback"]
    legend_row = length(channels) + 1
    CairoMakie.Legend(figure[legend_row, 1:3],
        vcat(size_handles, protocol_handles), vcat(size_labels, protocol_labels);
        nbanks=min(6, length(size_labels) + 3), framevisible=false,
        tellheight=true)
    policy = "Slice-local feedback: X uses one-step target feed-forward " *
        "(routing=$(scan.x_routing)); Z uses a complete three-sheet moving window. " *
        "Full-prefix feedback decodes the growing prepared volume."
    coincidence = "Coincident curves share the same coordinates; no visual offset is applied, " *
        "and every protocol remains identified by line style."
    CairoMakie.Label(figure[legend_row + 1, 1:3],
        "Individual logicals and union state failure · Paired differences use the same sampled faults\n" *
        policy * "\n" * coincidence;
        fontsize=12, tellwidth=false)
    return figure
end

function _toric3d_diagonal_feedback_save_figure(figure, path::AbstractString)
    extension = splitext(path)[2]
    temporary = path * ".tmp" * extension
    try
        CairoMakie.save(temporary, figure)
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

function _toric3d_diagonal_feedback_validate_for_output(scan, stem)
    invalid = scan.invalid_recoveries + sum(point.invalid_recoveries
        for point in scan.points)
    invalid == 0 && return nothing
    failure = stem * "-failure.toml"
    _toric3d_atomic_write(failure) do io
        TOML.print(io, Dict{String,Any}(
            "stage" => "artifact_validation",
            "invalid_recoveries" => scan.invalid_recoveries,
            "point_invalid_recoveries" => sum(point.invalid_recoveries
                for point in scan.points),
            "error" => "refusing to save a successful artifact with invalid recoveries");
            sorted=true)
    end
    error("refusing to save diagonal-feedback artifacts with invalid recoveries")
end

"""Save raw CSV, complete metadata, and PDF/SVG/PNG comparison figures."""
function save_toric3d_diagonal_feedback(scan::Toric3DDiagonalFeedbackScan,
        directory::AbstractString;
        basename::AbstractString=_TORIC3D_DIAGONAL_FEEDBACK_BASENAME)::NamedTuple
    stem = _toric3d_diagonal_feedback_output_stem(directory, basename)
    _toric3d_diagonal_feedback_validate_for_output(scan, stem)
    isempty(scan.points) && throw(ArgumentError("scan has no completed points"))
    csv = _write_toric3d_diagonal_feedback_raw(scan, stem * "-raw.csv")
    metadata = stem * "-metadata.toml"
    _toric3d_atomic_write(metadata) do io
        TOML.print(io, _toric3d_diagonal_feedback_metadata(scan); sorted=true)
    end
    figure = plot_toric3d_diagonal_feedback(scan)
    svg = _toric3d_diagonal_feedback_save_figure(figure, stem * ".svg")
    pdf = _toric3d_diagonal_feedback_save_figure(figure, stem * ".pdf")
    png = _toric3d_diagonal_feedback_save_figure(figure, stem * ".png")
    return (csv=csv, metadata=metadata, svg=svg, pdf=pdf, png=png)
end
