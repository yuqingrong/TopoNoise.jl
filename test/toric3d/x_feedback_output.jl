using Test
using Random
using TopoNoise

@testset "Diagonal X-feedback artifacts preserve paired provenance" begin
    scan = compare_toric3d_x_feedback_threshold(MersenneTwister(1234);
        sizes=[(2, 2, 2)], error_rates=[0.0, 0.01], shots=4,
        seed=1234, batch_size=3, routing=:one_step)

    mktempdir() do directory
        paths = save_toric3d_x_feedback_threshold(scan, directory)
        @test keys(paths) == (:csv, :metadata, :svg, :pdf, :png)
        @test all(path -> isfile(path) && filesize(path) > 0, values(paths))

        rows = readlines(paths.csv)
        @test length(rows) == 1 + 2 * 3
        header = split(first(rows), ',')
        required = [
            "steps", "Ly", "Lz", "qubits", "p", "metric", "shots",
            "seed", "point_seed", "routing", "final_only_failures",
            "slice_feedback_failures", "both_failures", "final_only_only",
            "slice_feedback_only", "final_only_rate", "final_only_ci_low",
            "final_only_ci_high", "slice_feedback_rate",
            "slice_feedback_ci_low", "slice_feedback_ci_high",
            "difference_feedback_minus_final", "difference_standard_error",
            "difference_ci_low", "difference_ci_high", "elapsed_seconds",
            "ldpc_version", "max_iter", "bp_method", "ms_scaling_factor",
            "schedule", "omp_thread_count", "osd_method", "osd_order"]
        @test all(in(header), required)
        metric_names = strip.(getindex.(split.(rows[2:end], ','), 6), '"')
        @test count(==("logical_1"), metric_names) == 2
        @test count(==("logical_2"), metric_names) == 2
        @test count(==("state_failure"), metric_names) == 2

        zero_rows = filter(row -> split(row, ',')[5] == "0.0", rows[2:end])
        @test length(zero_rows) == 3
        final_index = findfirst(==("final_only_failures"), header)
        feedback_index = findfirst(==("slice_feedback_failures"), header)
        @test all(row -> split(row, ',')[final_index] == "0" &&
            split(row, ',')[feedback_index] == "0", zero_rows)

        metadata = TopoNoise.TOML.parsefile(paths.metadata)
        @test contains(metadata["initial_state"], "|00>")
        @test metadata["logical_qubits"] == 2
        @test metadata["active_gate_outputs_per_transfer_cell"] == 13
        @test metadata["seed_gate_fault_locations"] == 0
        @test Set(metadata["ideal_operations"]) == Set([
            "initialization", "idle qubits", "syndrome measurements",
            "classical processing", "Pauli corrections"])
        @test contains(metadata["decoder_limitations"], "nonuniform exposure")
        @test contains(metadata["decoder_limitations"], "correlations")
        @test metadata["completed_points"] == 2
        @test metadata["expected_points"] == 2
        @test metadata["complete"]
        @test metadata["routing"] == "one_step"
        @test metadata["shots_per_point"] == 4
        @test metadata["batch_size"] == 3
        @test metadata["checkpoint_filename"] ==
            "toric3d-x-feedback-checkpoint.toml"
        @test metadata["resume_supported"]
        @test metadata["decoder_settings"]["max_iter"] == 50
        @test metadata["decoder_settings"]["bp_method"] == "minimum_sum"
        @test metadata["decoder_settings"]["osd_method"] == "OSD_CS"
        @test metadata["decoder_settings"]["osd_order"] == 2
    end

    empty_scan = Toric3DXFeedbackScan(scan.dimensions, scan.error_rates,
        scan.shots, scan.seed, scan.batch_size, scan.routing, scan.settings,
        Toric3DXFeedbackPoint[], scan.rng_type, scan.point_seeds,
        scan.started_at_unix, scan.resumed, scan.checkpoint_path)
    @test_throws ArgumentError plot_toric3d_x_feedback_threshold(empty_scan)
    mktempdir() do directory
        @test_throws ArgumentError save_toric3d_x_feedback_threshold(
            scan, directory; basename="")
        @test_throws ArgumentError save_toric3d_x_feedback_threshold(
            scan, directory; basename="nested/output")
    end
end
