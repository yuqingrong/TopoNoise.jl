using Test
using CairoMakie
using TopoNoise

function _feedback_strength_fixture(; invalid_recoveries=0)
    counts = (7, 0, 0, 0, 0, 0, 0, 3)
    metric = TopoNoise._toric3d_triple_metric(10, counts)
    metrics = (logical_1=metric, logical_2=metric, state_failure=metric)
    versions = (pymatching="2.4.0", ldpc="2.4.1")
    phase_seconds = (decoder_build=0.25, sampling=0.5, replay=1.25)
    x_point = Toric3DDiagonalFeedbackPoint((2, 2, 2), 0.01, :x_only, 10,
        UInt64(0x11), metrics, versions, phase_seconds, invalid_recoveries)
    z_point = Toric3DDiagonalFeedbackPoint((2, 2, 2), 0.01, :z_only, 10,
        UInt64(0x22), metrics, versions, phase_seconds, invalid_recoveries)
    return Toric3DDiagonalFeedbackScan([(2, 2, 2)], [0.01],
        [:x_only, :z_only], 10, 1234, 4, :one_step,
        Toric3DDecoderSettings(), [x_point, z_point], "Random.Xoshiro",
        UInt64[0x11, 0x22], 1.0, false, nothing, nothing,
        invalid_recoveries, (x=1, z=1))
end

function _feedback_figure_strings(figure)
    strings = String[]
    for block in figure.content
        if block isa CairoMakie.Axis
            append!(strings, string.((block.title[], block.xlabel[], block.ylabel[])))
        elseif block isa CairoMakie.Label
            push!(strings, string(block.text[]))
        end
        hasproperty(block, :scene) || continue
        for plot in block.scene.plots
            if hasproperty(plot, :text)
                value = plot.text[]
                append!(strings, string.(value isa AbstractArray ? value : (value,)))
            end
            if hasproperty(plot, :label)
                label = plot.label[]
                label isa AbstractString && !isempty(label) && push!(strings, label)
            end
        end
    end
    return strings
end

const _FEEDBACK_STRENGTH_COLUMNS = [
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

@testset "Diagonal feedback artifacts preserve three-way statistics" begin
    scan = _feedback_strength_fixture()
    figure = plot_toric3d_diagonal_feedback(scan)
    @test figure isa CairoMakie.Figure
    strings = join(_feedback_figure_strings(figure), '\n')
    @test contains(strings, "X-only")
    @test contains(strings, "Z-only")
    @test contains(strings, "No feedback (final decode only)")
    @test contains(strings, "Slice-local feedback")
    @test contains(strings, "Full-prefix feedback")
    @test contains(strings, "Individual logicals")
    @test contains(strings, "Paired differences")
    @test contains(lowercase(strings), "coincident")

    mktempdir() do directory
        paths = save_toric3d_diagonal_feedback(scan, directory)
        @test keys(paths) == (:csv, :metadata, :svg, :pdf, :png)
        @test all(path -> isfile(path) && filesize(path) > 100, values(paths))
        @test read(paths.pdf, 4) == UInt8[0x25, 0x50, 0x44, 0x46]
        @test read(paths.png, 8) == UInt8[0x89, 0x50, 0x4e, 0x47,
            0x0d, 0x0a, 0x1a, 0x0a]
        @test startswith(read(paths.svg, String), "<?xml")
        @test isempty(filter(name -> contains(name, ".tmp"), readdir(directory)))

        rows = readlines(paths.csv)
        @test length(rows) == 1 + 2 * 3
        header = split(first(rows), ',')
        @test header == _FEEDBACK_STRENGTH_COLUMNS
        indices = Dict(name => index for (index, name) in enumerate(header))
        keys_seen = Set{Tuple{String,String}}()
        for row in rows[2:end]
            fields = strip.(split(row, ','), '"')
            channel = fields[indices["channel"]]
            metric_name = fields[indices["metric"]]
            push!(keys_seen, (channel, metric_name))
            joint_range = indices["joint_n0_s0_f0"]:indices["joint_n1_s1_f1"]
            joint = parse.(Int, fields[joint_range])
            @test sum(joint) == 10
            @test parse(Int, fields[indices["none_failures"]]) ==
                sum(joint[mask + 1] for mask in 0:7 if (mask & 1) != 0)
            @test parse(Int, fields[indices["slice_local_failures"]]) ==
                sum(joint[mask + 1] for mask in 0:7 if (mask & 2) != 0)
            @test parse(Int, fields[indices["full_prefix_failures"]]) ==
                sum(joint[mask + 1] for mask in 0:7 if (mask & 4) != 0)
            @test parse(Float64, fields[indices[
                "none_vs_full_prefix_difference_full_prefix_minus_none"]]) == 0.0
            @test parse(Float64, fields[indices["decoder_build_seconds"]]) == 0.25
            @test fields[indices["pymatching_version"]] == "2.4.0"
        end
        @test keys_seen == Set((channel, metric) for channel in ("x_only", "z_only")
            for metric in ("logical_1", "logical_2", "state_failure"))

        metadata = TopoNoise.TOML.parsefile(paths.metadata)
        @test metadata["complete"]
        @test metadata["completed_points"] == metadata["expected_points"] == 2
        @test metadata["invalid_recoveries"] == 0
        @test Set(metadata["protocols"]) == Set(["none", "slice_local", "full_prefix"])
        @test metadata["joint_pattern_encoding"] ==
            "mask = none + 2*slice_local + 4*full_prefix"
        @test metadata["decoder_settings"]["max_iter"] == 50
        @test metadata["channel_models"]["x_only"]["seed_gate_outputs"] == 0
        @test metadata["channel_models"]["z_only"]["seed_gate_outputs_formula"] ==
            "3*Ly*Lz"
        x_locations = only(entry for entry in metadata["fault_locations"]
            if entry["channel"] == "x_only")
        z_locations = only(entry for entry in metadata["fault_locations"]
            if entry["channel"] == "z_only")
        @test x_locations["seed_gate_outputs"] == 0
        @test x_locations["transfer_gate_outputs"] == 104
        @test x_locations["total_gate_outputs"] == 104
        @test z_locations["seed_gate_outputs"] == 12
        @test z_locations["transfer_gate_outputs"] == 104
        @test z_locations["total_gate_outputs"] == 116
        @test metadata["source_hash_algorithm"] == "FNV-1a-64"
        @test !isempty(metadata["source_hashes"])
        @test all(occursin(r"^[0-9a-f]{16}$", value)
            for value in values(metadata["source_hashes"]))
        @test !isempty(metadata["git_state"]["commit"])
        @test haskey(metadata["git_state"], "dirty")
    end

    empty_scan = Toric3DDiagonalFeedbackScan(scan.dimensions, scan.error_rates,
        scan.channels, scan.shots, scan.seed, scan.batch_size, scan.x_routing,
        scan.settings, Toric3DDiagonalFeedbackPoint[], scan.rng_type,
        scan.point_seeds, scan.started_at_unix, scan.resumed,
        scan.checkpoint_path, scan.failure_path, scan.invalid_recoveries,
        scan.decoder_builds)
    @test_throws ArgumentError plot_toric3d_diagonal_feedback(empty_scan)

    mktempdir() do directory
        @test_throws ArgumentError save_toric3d_diagonal_feedback(
            scan, directory; basename="")
        @test_throws ArgumentError save_toric3d_diagonal_feedback(
            scan, directory; basename="nested/output")
    end

    invalid_scan = _feedback_strength_fixture(invalid_recoveries=1)
    mktempdir() do directory
        @test_throws ErrorException save_toric3d_diagonal_feedback(
            invalid_scan, directory)
        failure = joinpath(directory, "toric3d-diagonal-feedback-failure.toml")
        @test isfile(failure)
        payload = TopoNoise.TOML.parsefile(failure)
        @test payload["invalid_recoveries"] == 1
        @test payload["stage"] == "artifact_validation"
    end
end
