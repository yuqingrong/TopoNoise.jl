using Test

const _ISOMETRIC_PLANAR_CODE_EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "examples", "scan_isometric_planar_threshold.jl")

@testset "isometric planar-code threshold example" begin
    @test isfile(_ISOMETRIC_PLANAR_CODE_EXAMPLE_PATH)
    if isfile(_ISOMETRIC_PLANAR_CODE_EXAMPLE_PATH)
        include(_ISOMETRIC_PLANAR_CODE_EXAMPLE_PATH)
        defaults = ScanIsometricPlanarThresholdExample._parse_options(String[])
        @test defaults.distances == [3, 5, 7, 9]
        @test defaults.rates == collect(0.0:0.005:0.16)
        @test defaults.shots == 100_000
        @test defaults.batches == 100
        @test defaults.output_dir == normpath(joinpath(
            @__DIR__, "..", "results", "isometric-planar-code"))
        mktempdir() do directory
            output = IOBuffer()
            status = ScanIsometricPlanarThresholdExample.main([
                "--p-min", "0.08", "--p-max", "0.10", "--p-step", "0.02",
                "--distances", "3,5", "--shots", "20", "--batches", "2",
                "--bootstrap", "8", "--seed", "41", "--output-dir", directory];
                io=output)
            @test status == 0
            @test all(isfile(joinpath(directory, name)) for name in (
                "isometric_planar_capacity_scan.csv",
                "isometric_planar_capacity_batches.csv",
                "isometric_planar_capacity_crossings.csv",
                "isometric_planar_capacity_scaling_fit.csv",
                "isometric_planar_capacity_scaling_bootstrap.csv",
                "isometric_planar_capacity_scaling_loss_surface.csv",
                "isometric_planar_capacity_scaling_optimizer_trace.csv",
                "isometric_planar_capacity_scaling_sensitivity.csv",
                "isometric_planar_capacity_scan.svg",
                "isometric_planar_capacity_scan.pdf",
                "isometric_planar_capacity_scan.png",
                "isometric_planar_capacity_scaling_diagnostics.svg",
                "isometric_planar_capacity_scaling_diagnostics.pdf",
                "isometric_planar_capacity_scaling_diagnostics.png"))
            @test startswith(read(joinpath(
                directory, "isometric_planar_capacity_scan.svg"), String), "<?xml")
            @test startswith(read(joinpath(
                directory, "isometric_planar_capacity_scaling_diagnostics.svg"),
                String), "<?xml")
            @test startswith(read(joinpath(
                directory, "isometric_planar_capacity_scaling_diagnostics.pdf"),
                String), "%PDF")
            fit_csv = read(joinpath(
                directory, "isometric_planar_capacity_scaling_fit.csv"), String)
            @test startswith(fit_csv, "p_c,nu,p_c_se,nu_se,")
            @test occursin("insufficient_sizes", fit_csv)
            @test startswith(read(joinpath(directory,
                "isometric_planar_capacity_batches.csv"), String),
                "d,p,batch,batch_shots,failure_count,failure_rate")
            @test startswith(read(joinpath(directory,
                "isometric_planar_capacity_scaling_bootstrap.csv"), String),
                "replicate,status,input_rms_shift")
            @test startswith(read(joinpath(directory,
                "isometric_planar_capacity_scaling_loss_surface.csv"), String),
                "p_c,nu,loss,delta_loss,log10_one_plus_delta_loss")
            @test startswith(read(joinpath(directory,
                "isometric_planar_capacity_scaling_optimizer_trace.csv"), String),
                "step,p_c,nu,loss")
            @test startswith(read(joinpath(directory,
                "isometric_planar_capacity_scaling_sensitivity.csv"), String),
                "label,status,p_c,nu,loss")
        end
    end
end
