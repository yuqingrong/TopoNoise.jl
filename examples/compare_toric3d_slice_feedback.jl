module Toric3DSliceFeedbackExample

using Random
using TopoNoise

const USAGE = """
Usage: julia --project=. examples/compare_toric3d_slice_feedback.jl [options]

Compare unassisted transfer and slice matching feedback with identical gate faults.
Each size L means (steps,Ly,Lz)=(L,L,L). No final global recovery is applied.

Options:
  --sizes LIST        comma-separated cube sizes (default: 3,5,7)
  --error-rates LIST  comma-separated rates (default: 0,.001,.003,.01,.03,.05)
  --shots INTEGER     shots per channel point (default: 200)
  --seed INTEGER      nonnegative RNG seed (default: 1234)
  --routing NAME      one_step or full_path (default: one_step)
  --output-dir PATH   output directory (default: results/toric3d-slice-feedback-2026-10-11)
  --pilot-shots N     timing-pilot shots per channel (default: 5)
  --help              show this help
"""

function _list(value,parser,option)
    entries = split(value,',';keepempty=true)
    any(x -> isempty(strip(x)),entries) && throw(ArgumentError("$option requires a list without empty entries"))
    try
        return parser.(strip.(entries))
    catch
        throw(ArgumentError("$option has an invalid entry: $value"))
    end
end

function _parse_arguments(arguments)
    options = Dict{Symbol,Any}(:sizes=>[3,5,7],:error_rates=>[0.,.001,.003,.01,.03,.05],
        :shots=>200,:seed=>1234,:routing=>:one_step,:pilot_shots=>5,
        :output_dir=>joinpath("results","toric3d-slice-feedback-2026-10-11"))
    named = Dict("--sizes"=>:sizes,"--error-rates"=>:error_rates,"--shots"=>:shots,
        "--seed"=>:seed,"--routing"=>:routing,"--output-dir"=>:output_dir,"--pilot-shots"=>:pilot_shots)
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option=="--help" && return nothing
        haskey(named,option) || throw(ArgumentError("unknown option: $option"))
        index < length(arguments) && !startswith(arguments[index+1],"--") ||
            throw(ArgumentError("$option requires a value"))
        value,key = arguments[index+1],named[option]
        options[key] = if key===:sizes
            _list(value,x->parse(Int,x),option)
        elseif key===:error_rates
            _list(value,x->parse(Float64,x),option)
        elseif key===:output_dir
            String(value)
        elseif key===:routing
            Symbol(value)
        else
            parsed = tryparse(Int,value)
            parsed===nothing && throw(ArgumentError("$option requires an integer"))
            parsed
        end
        index += 2
    end
    sizes,rates = options[:sizes],options[:error_rates]
    all(x->x>=2,sizes) && all(diff(sizes).>0) || throw(ArgumentError("sizes must be strictly increasing integers at least 2"))
    all(x->isfinite(x) && 0<=x<.5,rates) && all(diff(rates).>0) || throw(ArgumentError("rates must be strictly increasing in [0,0.5)"))
    options[:shots]>0 && options[:pilot_shots]>0 || throw(ArgumentError("shot counts must be positive"))
    options[:seed]>=0 || throw(ArgumentError("seed must be nonnegative"))
    options[:routing] in (:one_step,:full_path) || throw(ArgumentError("routing must be one_step or full_path"))
    isempty(strip(options[:output_dir])) && throw(ArgumentError("output directory must not be empty"))
    return options
end

function main(arguments=ARGS;io::IO=stdout,error_io::IO=stderr)
    options = try
        _parse_arguments(arguments)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io,"error: ",sprint(showerror,error))
        println(error_io,USAGE)
        return 1
    end
    if options===nothing
        print(io,USAGE)
        return 0
    end
    seed,routing = options[:seed],options[:routing]
    try
        println(io,"Warmup: L=2, p=0.03, one paired shot per channel...")
        warmup_seconds = @elapsed compare_toric3d_slice_feedback(
            MersenneTwister(UInt(seed) ⊻ UInt(0x736c6963));sizes=[2],error_rates=[.03],shots=1,routing)
        println(io,"Warmup including compilation: ",round(warmup_seconds;digits=2)," s")
        pilot_size = maximum(options[:sizes])
        pilot_shots = options[:pilot_shots]
        println(io,"Timing pilot: L=$pilot_size, p=0.03, $pilot_shots paired shots per channel...")
        pilot_seconds = @elapsed compare_toric3d_slice_feedback(
            MersenneTwister(UInt(seed) ⊻ UInt(0x70696c6f));sizes=[pilot_size],error_rates=[.03],shots=pilot_shots,routing)
        combinations = length(options[:sizes])*length(options[:error_rates])
        estimate = pilot_seconds*options[:shots]*combinations/pilot_shots
        println(io,"Pilot runtime: ",round(pilot_seconds;digits=2)," s; rough scan estimate: ",
            round(estimate/60;digits=1)," min (size and matching cost vary).")
        println(io,"Running scan; completed points are saved in ",options[:output_dir])
        comparison = compare_toric3d_slice_feedback(MersenneTwister(seed);
            sizes=options[:sizes],error_rates=options[:error_rates],shots=options[:shots],
            seed,routing,output_dir=options[:output_dir],progress_io=io)
        paths = save_toric3d_slice_comparison(comparison,options[:output_dir])
        println(io,"Saved slice comparison: ",paths)
        println(io,"These physical syndrome diagnostics do not establish a logical-error threshold.")
    catch error
        println(error_io,"error: ",sprint(showerror,error))
        return 1
    end
    return 0
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(Toric3DSliceFeedbackExample.main())
end
