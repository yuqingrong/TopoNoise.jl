const _TRAJECTORY_EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "examples", "scan_toric_trajectories.jl")

@testset "Trajectory critical-scan CLI" begin
    @test isfile(_TRAJECTORY_EXAMPLE_PATH)

    if isfile(_TRAJECTORY_EXAMPLE_PATH)
        include(_TRAJECTORY_EXAMPLE_PATH)

        defaults = ScanToricTrajectoriesExample._parse_options([
            "--p-min", "0", "--p-max", "1", "--p-step", "1"])
        @test defaults.output_dir == normpath(joinpath(
            @__DIR__, "..", "results", "toric-trajectories"))

        mktempdir() do directory
            output = IOBuffer()
            arguments = [
                "--p-min", "0", "--p-max", "1", "--p-step", "1",
                "--sizes", "2,3", "--shots", "8", "--batches", "2",
                "--bootstrap", "4", "--seed", "1234",
                "--output-dir", directory,
            ]
            @test ScanToricTrajectoriesExample.main(
                arguments; io=output, error_io=output) == 0
            text = String(take!(output))

            scan_path = joinpath(directory, "trajectory_scan.csv")
            crossing_path = joinpath(directory, "trajectory_crossings.csv")
            for path in (scan_path, crossing_path)
                @test isfile(path)
                @test filesize(path) > 100
                @test contains(text, "wrote: $(abspath(path))")
            end
            for extension in ("svg", "pdf", "png")
                path = joinpath(directory, "trajectory_scan.$extension")
                @test isfile(path)
                @test filesize(path) > 100
                @test contains(text, "wrote: $(abspath(path))")
            end

            scan_lines = readlines(scan_path)
            @test length(scan_lines) == 5
            @test first(scan_lines) == join((
                "L", "p", "shots",
                "sampled_error_mean", "sampled_error_se",
                "boundary_error_mean", "boundary_error_se",
                "mismatch_mean", "mismatch_se",
                "frustration_mean", "frustration_se",
                "largest_cluster_mean", "largest_cluster_se",
                "horizontal_spanning_mean", "horizontal_spanning_se",
                "vertical_spanning_mean", "vertical_spanning_se"), ",")
            crossing_lines = readlines(crossing_path)
            @test length(crossing_lines) == 2
            @test first(crossing_lines) == join((
                "small_size", "large_size", "estimate", "ci_low",
                "ci_high", "confidence", "valid_bootstrap_fraction",
                "status"), ",")
            @test last(split(last(crossing_lines), ',')) in
                  ("ok", "unbracketed", "unstable")
        end

        for arguments in (
                String[],
                ["--p-min", "0", "--p-max", "1"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "0"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "0.3"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "1",
                 "--sizes", "4,2"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "1",
                 "--shots", "10", "--batches", "1"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "1",
                 "--spin-samples", "1"],
                ["--p-min", "0", "--p-max", "1", "--p-step", "1",
                 "--unknown", "value"])
            output = IOBuffer()
            @test ScanToricTrajectoriesExample.main(
                arguments; io=output, error_io=output) == 1
            text = String(take!(output))
            @test contains(text, "Usage:")
        end
    end
end
