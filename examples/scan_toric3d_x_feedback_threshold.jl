module Toric3DXFeedbackExample

using Random
using TopoNoise
import TOML

const USAGE = """
Usage: julia --project=. examples/scan_toric3d_x_feedback_threshold.jl [options]
  --sizes LIST                    cubic sizes, default 3,5,7,9
  --error-rates LIST              default 0,.0001,.0003,.001,.002,.003,.005,.0075,.01,.015,.02
  --shots N                       default 1000 paired shots per point
  --seed N                        default 1234
  --batch-size N                  default 256 checkpoint batch size
  --pilot-shots N                 default 20 independent timing shots
  --routing one_step|full_path    default one_step
  --output-dir PATH               default results/toric3d-x-feedback-threshold-2026-10-11/pilot
  --resume                        resume the matching checkpoint
  --help

X-only post-H/CNOT output faults in the encoded diagonal sweep. The seed is
an ideally prepared logical |00> state. Measurements, processing, corrections,
initialization, and idle qubits are ideal. Failure means either residual
logical-X membrane after validated final BP+OSD recovery.
"""

const _DEFAULTS = Dict(
    "--sizes" => "3,5,7,9",
    "--error-rates" => "0,.0001,.0003,.001,.002,.003,.005,.0075,.01,.015,.02",
    "--shots" => "1000",
    "--seed" => "1234",
    "--batch-size" => "256",
    "--pilot-shots" => "20",
    "--routing" => "one_step",
    "--output-dir" => "results/toric3d-x-feedback-threshold-2026-10-11/pilot")

function _parse_list(::Type{T}, text::AbstractString, name::AbstractString) where T
    isempty(strip(text)) && throw(ArgumentError("$name must not be empty"))
    values = try
        parse.(T, split(text, ','))
    catch error
        throw(ArgumentError("could not parse $name: $(sprint(showerror, error))"))
    end
    return values
end

function _parse_arguments(arguments)
    options = copy(_DEFAULTS)
    resume = false
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option == "--resume" && begin
            resume = true
            index += 1
            continue
        end
        haskey(options, option) || throw(ArgumentError("unknown option $option"))
        index < length(arguments) || throw(ArgumentError("$option needs a value"))
        value = arguments[index + 1]
        startswith(value, "--") && throw(ArgumentError("$option needs a value"))
        options[option] = value
        index += 2
    end

    sizes = _parse_list(Int, options["--sizes"], "sizes")
    rates = _parse_list(Float64, options["--error-rates"], "error rates")
    shots = parse(Int, options["--shots"])
    seed = parse(Int, options["--seed"])
    batch_size = parse(Int, options["--batch-size"])
    pilot_shots = parse(Int, options["--pilot-shots"])
    all(size -> size >= 2, sizes) || throw(ArgumentError("sizes must be at least 2"))
    issorted(sizes; lt=<) || throw(ArgumentError("sizes must be strictly increasing"))
    all(rate -> isfinite(rate) && 0 <= rate < 0.5, rates) ||
        throw(ArgumentError("error rates must lie in [0, 0.5)"))
    issorted(rates; lt=<) || throw(ArgumentError("error rates must be strictly increasing"))
    shots > 0 || throw(ArgumentError("shots must be positive"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    batch_size > 0 || throw(ArgumentError("batch size must be positive"))
    pilot_shots > 0 || throw(ArgumentError("pilot shots must be positive"))
    routing = Symbol(options["--routing"])
    routing in (:one_step, :full_path) ||
        throw(ArgumentError("routing must be one_step or full_path"))
    directory = options["--output-dir"]
    isempty(strip(directory)) && throw(ArgumentError("output directory must not be empty"))
    return (; sizes, error_rates=rates, shots, seed, batch_size, pilot_shots,
        routing, output_dir=directory, resume)
end

function _write_timing(path, data)
    temporary = path * ".tmp"
    open(temporary, "w") do timing_io
        TOML.print(timing_io, data; sorted=true)
    end
    mv(temporary, path; force=true)
    return path
end

function _timing_pilot(options; io::IO)
    warmup_rng = MersenneTwister(UInt(options.seed) ⊻ UInt(0x7761726d))
    println(io, "Warmup: L=2, p=.01, two paired shots...")
    flush(io)
    warmup_seconds = @elapsed compare_toric3d_x_feedback_threshold(warmup_rng;
        sizes=[2], error_rates=[0.01], shots=2, batch_size=2,
        routing=options.routing)

    pilot_size = maximum(options.sizes)
    pilot_p = 0.01
    println(io, "Timing pilot: L=$pilot_size, p=$pilot_p, $(options.pilot_shots) shots...")
    flush(io)
    setup_seconds = @elapsed begin
        code = Toric3DEncodedSlice(pilot_size, pilot_size, pilot_size)
        encoder = toric3d_encoded_encoder(code)
        slice_decoder = build_toric3d_slice_decoder(code.lattice)
        final_decoder = build_toric3d_x_threshold_decoder(code; p=pilot_p)
    end
    pilot_rng = MersenneTwister(UInt(options.seed) ⊻ UInt(0x70696c6f74))
    first_decode_seconds = @elapsed begin
        faults = sample_toric3d_faults(pilot_rng, encoder; p_x=pilot_p, p_z=0.0)
        toric3d_x_feedback_pair(encoder, faults, slice_decoder, final_decoder;
            routing=options.routing)
    end
    pilot_seconds = @elapsed for _ in 1:options.pilot_shots
        faults = sample_toric3d_faults(pilot_rng, encoder; p_x=pilot_p, p_z=0.0)
        toric3d_x_feedback_pair(encoder, faults, slice_decoder, final_decoder;
            routing=options.routing)
    end
    point_count = length(options.sizes) * length(options.error_rates)
    estimated_scan_seconds = (setup_seconds + first_decode_seconds) * point_count +
        pilot_seconds / options.pilot_shots * options.shots * point_count
    println(io, "Warmup ", round(warmup_seconds; digits=2), " s; setup ",
        round(setup_seconds; digits=3), " s; first decode ",
        round(first_decode_seconds; digits=3), " s; timing shots ",
        round(pilot_seconds; digits=3), " s.")
    println(io, "Rough scan estimate: ", round(estimated_scan_seconds / 60; digits=1),
        " min.")
    flush(io)
    return Dict{String,Any}(
        "warmup_seconds" => warmup_seconds,
        "warmup_shots" => 2,
        "setup_seconds" => setup_seconds,
        "first_decode_seconds" => first_decode_seconds,
        "pilot_seconds" => pilot_seconds,
        "pilot_size" => pilot_size,
        "pilot_p" => pilot_p,
        "pilot_shots" => options.pilot_shots,
        "production_shots_per_point" => options.shots,
        "production_points" => point_count,
        "estimated_scan_seconds" => estimated_scan_seconds,
        "timing_shots_included_in_production" => false,
        "routing" => string(options.routing))
end

function main(arguments=ARGS; io::IO=stdout, error_io::IO=stderr)::Int
    if "--help" in arguments
        print(io, USAGE)
        return 0
    end
    try
        options = _parse_arguments(arguments)
        mkpath(options.output_dir)
        timing = _timing_pilot(options; io)
        _write_timing(joinpath(options.output_dir, "timing-pilot.toml"), timing)
        options.resume && println(io, "Resuming matching checkpoint in $(options.output_dir)...")
        scan = compare_toric3d_x_feedback_threshold(MersenneTwister(options.seed);
            sizes=options.sizes, error_rates=options.error_rates,
            shots=options.shots, seed=options.seed,
            batch_size=options.batch_size, routing=options.routing,
            output_dir=options.output_dir, resume=options.resume, progress_io=io)
        paths = save_toric3d_x_feedback_threshold(scan, options.output_dir)
        println(io, "Saved: ", paths)
        flush(io)
        return 0
    catch error
        error isa InterruptException && rethrow()
        print(error_io, "error: ")
        showerror(error_io, error)
        println(error_io)
        print(error_io, USAGE)
        return 1
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(Toric3DXFeedbackExample.main())
end
