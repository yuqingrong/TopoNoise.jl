module Toric3DZAspectRatioExample

using Random
using TopoNoise
import TOML

const USAGE = """
Usage: julia --project=. examples/compare_toric3d_z_aspect_ratio.jl [options]
  --sizes LIST                    transverse sizes L, default 3,5,7
  --step-multipliers LIST         step counts T/L, default 1,2
  --error-rates LIST              default .0045,.005,.0055,.006,.0065
  --shots N                       default 1000 per point
  --seed N                        default 1234
  --batch-size N                  default 256
  --pilot-shots N                 default 20 independent timing shots
  --output-dir PATH               default results/toric3d-h6-z-aspect-ratio-2026-10-15/production
  --help

Z-only output faults in the plaquette-transfer H+6-CNOT circuit. This is a
T=L versus T=2L aspect-ratio diagnostic with final-only decoding and no intermediate feedback.
"""

const _DEFAULTS = Dict(
    "--sizes" => "3,5,7",
    "--step-multipliers" => "1,2",
    "--error-rates" => ".0045,.005,.0055,.006,.0065",
    "--shots" => "1000",
    "--seed" => "1234",
    "--batch-size" => "256",
    "--pilot-shots" => "20",
    "--output-dir" =>
        "results/toric3d-h6-z-aspect-ratio-2026-10-15/production")

function _parse_list(::Type{T}, text::AbstractString, name::AbstractString) where T
    isempty(strip(text)) && throw(ArgumentError("$name must not be empty"))
    return try
        parse.(T, split(text,','))
    catch error
        throw(ArgumentError("could not parse $name: $(sprint(showerror,error))"))
    end
end

function _positive_unique_integers(values, name::AbstractString)
    all(value -> value isa Integer && value > 0,values) ||
        throw(ArgumentError("$name must contain only positive integers"))
    length(unique(values)) == length(values) ||
        throw(ArgumentError("$name must be unique"))
    return values
end

function _aspect_dimensions(sizes,step_multipliers)
    _positive_unique_integers(sizes,"sizes")
    _positive_unique_integers(step_multipliers,"step multipliers")
    dimensions = NTuple{3,Int}[]
    for L in sizes, multiplier in step_multipliers
        steps = try
            Base.checked_mul(L,multiplier)
        catch error
            error isa OverflowError || rethrow()
            throw(ArgumentError("size times step multiplier overflows Int"))
        end
        push!(dimensions,(steps,L,L))
    end
    return TopoNoise._toric3d_sizes(dimensions)
end

function _parse_arguments(arguments)
    options = copy(_DEFAULTS)
    iseven(length(arguments)) || throw(ArgumentError("every option needs a value"))
    for index in 1:2:length(arguments)
        option = arguments[index]
        haskey(options,option) || throw(ArgumentError("unknown option $option"))
        value = arguments[index+1]
        startswith(value,"--") && throw(ArgumentError("$option needs a value"))
        options[option] = value
    end
    sizes = _parse_list(Int,options["--sizes"],"sizes")
    step_multipliers = _parse_list(Int,options["--step-multipliers"],"step multipliers")
    dimensions = _aspect_dimensions(sizes,step_multipliers)
    rates = _parse_list(Float64,options["--error-rates"],"error rates")
    error_rates = TopoNoise._scan_error_rates(rates)
    shots = parse(Int,options["--shots"])
    seed = parse(Int,options["--seed"])
    batch_size = parse(Int,options["--batch-size"])
    pilot_shots = parse(Int,options["--pilot-shots"])
    shots > 0 || throw(ArgumentError("shots must be positive"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    batch_size > 0 || throw(ArgumentError("batch size must be positive"))
    pilot_shots > 0 || throw(ArgumentError("pilot shots must be positive"))
    output_dir = options["--output-dir"]
    isempty(strip(output_dir)) && throw(ArgumentError("output directory must not be empty"))
    return (;sizes,step_multipliers,dimensions,error_rates,shots,seed,batch_size,
        pilot_shots,output_dir)
end

_qubit_count((steps,Ly,Lz)) = 3(steps+1)*Ly*Lz
_fault_location_count((steps,Ly,Lz)) = (3+13steps)*Ly*Lz

function _timing_target(dimensions)
    checked = TopoNoise._toric3d_sizes(dimensions)
    _,index = findmax(dimension ->
        (_fault_location_count(dimension),_qubit_count(dimension)),checked)
    return checked[index]
end

function _atomic_toml(path::AbstractString,data)
    temporary = path*".tmp"
    open(temporary,"w") do output
        TOML.print(output,data;sorted=true)
    end
    mv(temporary,path;force=true)
    return path
end

function _timing_pilot(options;io::IO)
    println(io,"Warmup: T=1, Ly=Lz=2, p=.01, two shots...")
    flush(io)
    warmup_seconds = @elapsed scan_toric3d_z_threshold(
        MersenneTwister(xor(UInt(options.seed),UInt(0x7761726d)));
        sizes=[(1,2,2)],error_rates=[.01],shots=2,batch_size=2)

    pilot_dimensions = _timing_target(options.dimensions)
    pilot_p = .01
    println(io,"Timing pilot: dimensions=$pilot_dimensions, p=$pilot_p, ",
        options.pilot_shots," shots...")
    flush(io)
    setup_seconds = @elapsed begin
        code = Toric3DEncodedSlice(pilot_dimensions...)
        model = build_toric3d_z_fault_model(code)
        decoder = build_toric3d_z_threshold_decoder(code)
    end
    pilot_rng = MersenneTwister(xor(UInt(options.seed),UInt(0x70696c6f74)))
    first_decode_seconds = @elapsed toric3d_z_threshold_counts(decoder,
        sample_toric3d_z_frames(pilot_rng,model;p=pilot_p,shots=1))
    pilot_seconds = @elapsed toric3d_z_threshold_counts(decoder,
        sample_toric3d_z_frames(pilot_rng,model;p=pilot_p,
            shots=options.pilot_shots))
    point_count = length(options.dimensions)*length(options.error_rates)
    estimated_scan_seconds = (setup_seconds+first_decode_seconds)*length(options.dimensions)+
        pilot_seconds/options.pilot_shots*options.shots*point_count
    println(io,"Warmup ",round(warmup_seconds;digits=2)," s; setup ",
        round(setup_seconds;digits=3)," s; first decode ",
        round(first_decode_seconds;digits=3)," s; timing shots ",
        round(pilot_seconds;digits=3)," s.")
    println(io,"Rough scan estimate: ",round(estimated_scan_seconds/60;digits=1)," min.")
    flush(io)
    return Dict{String,Any}(
        "warmup_seconds"=>warmup_seconds,
        "warmup_shots"=>2,
        "setup_seconds"=>setup_seconds,
        "first_decode_seconds"=>first_decode_seconds,
        "pilot_seconds"=>pilot_seconds,
        "pilot_dimensions_steps_Ly_Lz"=>collect(pilot_dimensions),
        "pilot_qubits"=>_qubit_count(pilot_dimensions),
        "pilot_fault_locations"=>_fault_location_count(pilot_dimensions),
        "pilot_p"=>pilot_p,
        "pilot_shots"=>options.pilot_shots,
        "production_shots_per_point"=>options.shots,
        "production_points"=>point_count,
        "estimated_scan_seconds"=>estimated_scan_seconds,
        "timing_shots_included_in_production"=>false)
end

function _study_metadata(options)
    return Dict{String,Any}(
        "format_version"=>1,
        "study"=>"H+6 Z-only aspect-ratio comparison",
        "circuit"=>"plaquette-transfer H+6-CNOT",
        "channel"=>"Z-only post-H/CNOT output faults",
        "decoding"=>"final-only",
        "intermediate_feedback"=>"none",
        "comparison"=>"T=L versus T=2L aspect-ratio diagnostic",
        "threshold_eligible"=>false,
        "dimensions_steps_Ly_Lz"=>[collect(dimension) for dimension in options.dimensions],
        "transverse_sizes"=>options.sizes,
        "step_multipliers"=>options.step_multipliers,
        "error_rates"=>options.error_rates,
        "shots_per_point"=>options.shots,
        "seed"=>options.seed,
        "batch_size"=>options.batch_size)
end

function main(arguments=ARGS;io::IO=stdout,error_io::IO=stderr)::Int
    if "--help" in arguments
        print(io,USAGE)
        return 0
    end
    try
        options = _parse_arguments(arguments)
        mkpath(options.output_dir)
        timing = _timing_pilot(options;io)
        _atomic_toml(joinpath(options.output_dir,"timing-pilot.toml"),timing)
        scan = scan_toric3d_z_threshold(MersenneTwister(options.seed);
            sizes=options.dimensions,error_rates=options.error_rates,
            shots=options.shots,seed=options.seed,batch_size=options.batch_size,
            output_dir=options.output_dir,progress_io=io)
        paths = save_toric3d_z_threshold(scan,options.output_dir)
        study = _atomic_toml(joinpath(options.output_dir,
            "toric3d-z-aspect-ratio-study.toml"),_study_metadata(options))
        println(io,"Saved: ",(core=paths,study=study))
        flush(io)
        return 0
    catch error
        error isa InterruptException && rethrow()
        print(error_io,"error: ")
        showerror(error_io,error)
        println(error_io)
        print(error_io,USAGE)
        return 1
    end
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(Toric3DZAspectRatioExample.main())
end
