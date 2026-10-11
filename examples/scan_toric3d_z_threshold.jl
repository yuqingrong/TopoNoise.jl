module Toric3DZThresholdExample
using Random, TopoNoise
import TOML

const USAGE = """
Usage: julia --project=. examples/scan_toric3d_z_threshold.jl [options]
  --sizes LIST         cubic dimensions, default 5,7,9
  --error-rates LIST   default .002,.004,.006,.008,.010,.015,.020
  --shots N           default 1000 per point
  --seed N            default 1234
  --batch-size N      default 256
  --pilot-shots N     default 20 for the timing pilot
  --output-dir PATH   default results/toric3d-z-threshold-2026-10-11/pilot
  --help

Z-only post-H/CNOT noise; noisy-H measurement-prepared logical |++> seed;
final unit-weight data-edge MWPM; failure is either residual logical Z.
"""

function main(arguments=ARGS)
    "--help" in arguments && (print(USAGE); return 0)
    options = Dict("--sizes"=>"5,7,9","--error-rates"=>".002,.004,.006,.008,.010,.015,.020",
        "--shots"=>"1000","--seed"=>"1234","--batch-size"=>"256","--pilot-shots"=>"20",
        "--output-dir"=>"results/toric3d-z-threshold-2026-10-11/pilot")
    try
        iseven(length(arguments)) || throw(ArgumentError("every option needs a value"))
        for index in 1:2:length(arguments)
            haskey(options,arguments[index]) || throw(ArgumentError("unknown option $(arguments[index])"))
            options[arguments[index]] = arguments[index+1]
        end
        sizes = parse.(Int,split(options["--sizes"],','))
        rates = parse.(Float64,split(options["--error-rates"],','))
        shots,seed,batch,pilot_shots = (parse(Int,options[key])
            for key in ("--shots","--seed","--batch-size","--pilot-shots"))
        seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
        directory = options["--output-dir"]
        println("Warmup: L=2, p=.01, two shots..."); flush(stdout)
        warmup = @elapsed scan_toric3d_z_threshold(MersenneTwister(UInt(seed) ⊻ UInt(0x7761726d));
            sizes=[2],error_rates=[.01],shots=2,batch_size=batch)
        println("Timing pilot: L=",maximum(sizes),", p=.01, ",pilot_shots," shots..."); flush(stdout)
        pilot_code = Toric3DEncodedSlice(maximum(sizes),maximum(sizes),maximum(sizes))
        setup = @elapsed begin
            pilot_model = build_toric3d_z_fault_model(pilot_code)
            pilot_decoder = build_toric3d_z_threshold_decoder(pilot_code)
        end
        # PyMatching builds its internal decoding structures on first use.
        # Account for that once per size, rather than multiplying it per shot.
        decoder_warmup = @elapsed toric3d_z_threshold_counts(pilot_decoder,
            sample_toric3d_z_frames(MersenneTwister(UInt(seed) ⊻ UInt(0x6d617463)),
                pilot_model;p=.01,shots=1))
        pilot = @elapsed toric3d_z_threshold_counts(pilot_decoder,
            sample_toric3d_z_frames(MersenneTwister(UInt(seed) ⊻ UInt(0x70696c6f)),
                pilot_model;p=.01,shots=pilot_shots))
        estimate = (setup+decoder_warmup)*length(sizes) +
            pilot*shots*length(sizes)*length(rates)/pilot_shots
        println("Warmup ",round(warmup;digits=2)," s; pilot ",round(pilot;digits=3),
            " s; rough scan estimate ",round(estimate/60;digits=1)," min."); flush(stdout)
        mkpath(directory)
        open(joinpath(directory,"timing-pilot.toml"),"w") do io
            TOML.print(io,Dict("warmup_seconds"=>warmup,"pilot_setup_seconds"=>setup,"pilot_seconds"=>pilot,
                "decoder_warmup_seconds"=>decoder_warmup,
                "pilot_size"=>maximum(sizes),"pilot_shots"=>pilot_shots,"pilot_p"=>.01,
                "estimated_scan_seconds"=>estimate);sorted=true)
        end
        scan = scan_toric3d_z_threshold(MersenneTwister(seed);sizes,error_rates=rates,
            shots,seed,batch_size=batch,output_dir=directory,progress_io=stdout)
        println("Saved: ",save_toric3d_z_threshold(scan,directory))
        return 0
    catch error
        showerror(stderr,error); println(stderr); return 1
    end
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(Toric3DZThresholdExample.main())
end
