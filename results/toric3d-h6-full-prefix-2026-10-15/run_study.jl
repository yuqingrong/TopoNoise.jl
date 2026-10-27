module H6FullPrefixStudy

using TopoNoise, Random, LinearAlgebra
import TOML, SHA, PythonCall

const ROOT=abspath(joinpath(@__DIR__,"..",".."))
const ARMS=("final_only","full_prefix")
const USAGE="""
Usage: julia --project=. --startup-file=no run_study.jl --channel x_only|z_only --output DIR --sizes 3,7,11 --rates 0,.005,.01 --shots N --seed N [--batch-size 25] [--resume] [--max-batches N]

Paired legacy H+6 final-only versus complete full-prefix physical recovery, T=L.
Feedback occurs after EVERY transfer 1:T, including T, before one final outcome.
Each channel uses one immutable fault record per pair and the same final decoder.
X starts from ideal |00>; Z has noisy seed-H |++> preparation without r=0 repair.
--max-batches stops after committing N batches; --resume restores the exact
point substreams and rebuilds raw CSV without appending duplicate rows.
"""

digest(bytes)=bytes2hex(SHA.sha256(bytes))
hexseed(seed)=string(seed;base=16,pad=16)

function parse_arguments(arguments)
    values=Dict("--batch-size"=>"25","--max-batches"=>string(typemax(Int)))
    required=("--channel","--output","--sizes","--rates","--shots","--seed")
    allowed=(required...,"--batch-size","--max-batches")
    seen=Set{String}(); resume=false; i=1
    while i <= length(arguments)
        key=arguments[i]
        key in seen && throw(ArgumentError("repeated option: $key"))
        push!(seen,key)
        if key == "--resume"
            resume=true; i+=1; continue
        end
        key in allowed || throw(ArgumentError("unknown option: $key"))
        i < length(arguments) || throw(ArgumentError("$key requires a value"))
        values[key]=arguments[i+1]; i+=2
    end
    all(haskey(values,key) for key in required) || throw(ArgumentError("missing required option"))
    channel=values["--channel"]
    channel in ("x_only","z_only") || throw(ArgumentError("channel must be x_only or z_only"))
    sizes=parse.(Int,split(values["--sizes"],','))
    rates=parse.(Float64,split(values["--rates"],','))
    !isempty(sizes) && all(>=(2),sizes) && all(i -> sizes[i]<sizes[i+1],1:length(sizes)-1) ||
        throw(ArgumentError("sizes must be strictly increasing integers >= 2"))
    !isempty(rates) && all(p -> isfinite(p) && 0 <= p < .5,rates) &&
        all(i -> rates[i]<rates[i+1],1:length(rates)-1) ||
        throw(ArgumentError("rates must be strictly increasing values in [0,.5)"))
    shots=parse(Int,values["--shots"]); seed=parse(Int,values["--seed"])
    batch_size=parse(Int,values["--batch-size"]); max_batches=parse(Int,values["--max-batches"])
    shots > 0 && seed >= 0 && batch_size > 0 && max_batches > 0 ||
        throw(ArgumentError("shots, batch-size, max-batches must be positive; seed must be nonnegative"))
    output=values["--output"]
    !isempty(strip(output)) || throw(ArgumentError("output directory must not be empty"))
    config=Dict{String,Any}("channel"=>channel,"sizes"=>sizes,"rates"=>rates,"shots"=>shots,
        "seed"=>seed,"batch_size"=>batch_size,"dimensions"=>[[L,L,L] for L in sizes])
    return (;config,output=abspath(output),resume,max_batches)
end

"""Size-major, rate-minor seeds: scalar UInt64 draws; duplicate values discarded."""
function point_seeds(seed::Integer,count::Integer)
    seed >= 0 && count > 0 || throw(ArgumentError("seed >= 0 and count > 0 required"))
    rng=MersenneTwister(seed); result=UInt64[]; seen=Set{UInt64}()
    while length(result) < count
        value=rand(rng,UInt64)
        value in seen && continue
        push!(seen,value); push!(result,value)
    end
    return result
end

failure_mask(flags)=sum(Int(flags[i]) << (i-1) for i in 1:2)
point_grid(config)=[(channel=config["channel"],L=L,p=p) for L in config["sizes"] for p in config["rates"]]

function sample_faults(rng,encoder,p,channel)
    channel == "x_only" && return sample_toric3d_faults(rng,encoder;p_x=p,p_z=0.0)
    channel == "z_only" && return sample_toric3d_diagonal_z_faults(rng,encoder;p)
    throw(ArgumentError("unsupported noise channel"))
end

function resume_rng(seed::UInt64,completed::Integer,encoder,p,channel)
    completed >= 0 || throw(ArgumentError("completed records must be nonnegative"))
    rng=Random.Xoshiro(seed)
    for _ in 1:completed
        sample_faults(rng,encoder,p,channel)
    end
    return rng
end

function toml_text(data)
    io=IOBuffer(); TOML.print(io,data;sorted=true); return String(take!(io))
end

function atomic_write(path,contents)
    temporary=path*".tmp"
    open(temporary,"w") do io
        write(io,contents)
    end
    mv(temporary,path;force=true)
    return path
end

function source_hashes()
    paths=["Project.toml",relpath(@__FILE__,ROOT)]
    isfile(joinpath(ROOT,"Manifest.toml")) && push!(paths,"Manifest.toml")
    for (directory,_,files) in walkdir(joinpath(ROOT,"src")), file in files
        endswith(file,".jl") && push!(paths,relpath(joinpath(directory,file),ROOT))
    end
    return Dict(path=>digest(read(joinpath(ROOT,path))) for path in sort!(paths))
end

function runtime_metadata()
    Threads.nthreads() == 1 && BLAS.get_num_threads() == 1 ||
        throw(ArgumentError("one Julia and BLAS thread required"))
    thread_keys=("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    all(get(ENV,key,"") == "1" for key in thread_keys) ||
        throw(ArgumentError("all four thread environment variables must equal 1"))
    get(ENV,"JULIA_CONDAPKG_BACKEND","") == "Null" || throw(ArgumentError("JULIA_CONDAPKG_BACKEND=Null required"))
    executable=PythonCall.pyconvert(String,PythonCall.pyimport("sys").executable)
    requested=get(ENV,"JULIA_PYTHONCALL_EXE","")
    !isempty(requested) && realpath(executable) == realpath(requested) ||
        throw(ArgumentError("Python runtime must equal JULIA_PYTHONCALL_EXE"))
    versions=PythonCall.pyimport("importlib.metadata")
    return Dict{String,Any}("julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),"python_executable"=>executable,
        "python_version"=>PythonCall.pyconvert(String,PythonCall.pyimport("sys").version),
        "packages"=>Dict(name=>PythonCall.pyconvert(String,versions.version(name)) for name in ("numpy","scipy","ldpc","pymatching")),
        "environment"=>Dict(key=>get(ENV,key,"") for key in (thread_keys...,"JULIA_CONDAPKG_BACKEND","JULIA_PYTHONCALL_EXE")))
end

function protocol_metadata(channel)
    settings=TopoNoise._toric3d_settings_metadata(Toric3DDecoderSettings())
    return Dict{String,Any}("circuit"=>"legacy H+6-CNOT plaquette transfer","channel"=>channel,
        "geometry"=>"T=L; Ly=Lz=L","arm_order"=>collect(ARMS),
        "initial_state"=>channel == "x_only" ? "ideal encoded logical |00>" : "noisy seed-H preparation of logical |++>",
        "seed_preparation"=>channel == "x_only" ? "ideal; zero seed fault locations" : "3 L^2 noisy seed-H outputs; no separate r=0 seed correction in either arm",
        "final_only"=>"no intermediate physical correction; exactly one final recovery and logical outcome",
        "full_prefix"=>"complete physical recovery on all available data slices Q0..Qr after EVERY transfer 1:T, including the last; complete-prefix syndrome zero; exactly one final recovery and logical outcome afterward",
        "feedback_timing"=>"after transfers 1:T inclusive; no separate r=0 correction",
        "fault_locations"=>channel == "x_only" ? "13 L^3 participating transfer outputs" : "3 L^2 seed-H outputs plus 13 L^3 participating transfer outputs",
        "noise"=>"independent Bernoulli-p faults in the selected channel at its recorded circuit locations",
        "sampler"=>channel == "x_only" ? "sample_toric3d_faults(p_x=p,p_z=0)" : "sample_toric3d_diagonal_z_faults: seed first, then participating transfer outputs in schedule order",
        "shared_final_decoder"=>channel == "x_only" ? "family.prefixes[end]; full-plaquette BP+OSD with declared p prior; syndrome-validated recovery" : "family.prefixes[end]; unit-weight data-edge PyMatching MWPM; syndrome-validated recovery",
        "decoder_settings"=>channel == "x_only" ? Dict(string(k)=>v for (k,v) in pairs(settings)) :
            Dict("weights"=>1.0,"merge_strategy"=>"smallest-weight"),
        "z_weights"=>"unit data-edge weights", "failure"=>"either residual logical parity after final recovery",
        "joint_histogram_encoding"=>"mask = final_only + 2*full_prefix; each term is its state-failure Boolean",
        "point_order"=>"sizes outer loop, rates inner loop; both strictly increasing",
        "point_seed_generation"=>"MersenneTwister(integer base seed); repeated scalar rand(rng,UInt64); discard duplicate values; 16 lowercase hexadecimal digits",
        "point_rng"=>"Random.Xoshiro(point_seed)",
        "resume"=>"recreate point Xoshiro and sample/discard completed records without decoding; rebuild raw CSV from authoritative checkpoint")
end

function fresh_checkpoint(config,runtime,hashes)
    count=length(point_grid(config))
    states=[Dict{String,Any}("completed_shots"=>0,"logical_1_counts"=>zeros(Int,2),
        "logical_2_counts"=>zeros(Int,2),"failures"=>zeros(Int,2),"joint_counts"=>zeros(Int,4),
        "elapsed_seconds"=>0.0) for _ in 1:count]
    return Dict{String,Any}("format_version"=>1,"config"=>config,"runtime"=>runtime,
        "source_hashes"=>hashes,"protocol"=>protocol_metadata(config["channel"]),"point_seeds_hex"=>hexseed.(point_seeds(config["seed"],count)),
        "states"=>states,"generation"=>0,"started_at_unix"=>time(),"updated_at_unix"=>time(),
        "decoder_build_seconds_by_size"=>Dict(string(L)=>0.0 for L in config["sizes"]),
        "decoder_build_count_by_size"=>Dict(string(L)=>0 for L in config["sizes"]))
end

function validate_states(checkpoint)
    config=checkpoint["config"]; states=checkpoint["states"]
    length(states) == length(point_grid(config)) || throw(ArgumentError("checkpoint grid length mismatch"))
    checkpoint["point_seeds_hex"] == hexseed.(point_seeds(config["seed"],length(states))) ||
        throw(ArgumentError("checkpoint point seeds differ from exact base-seed regeneration"))
    for L in config["sizes"]
        seconds=checkpoint["decoder_build_seconds_by_size"][string(L)]
        builds=checkpoint["decoder_build_count_by_size"][string(L)]
        isfinite(seconds) && seconds >= 0 && builds isa Integer && builds >= 0 ||
            throw(ArgumentError("invalid decoder build timing"))
    end
    for state in states
        n=state["completed_shots"]
        n isa Integer && !(n isa Bool) && 0 <= n <= config["shots"] || throw(ArgumentError("invalid completed sample count"))
        for (key,size) in (("logical_1_counts",2),("logical_2_counts",2),("failures",2),("joint_counts",4))
            counts=state[key]
            length(counts) == size && all(x -> x isa Integer && !(x isa Bool) && 0 <= x <= n,counts) ||
                throw(ArgumentError("invalid $key"))
        end
        sum(state["joint_counts"]) == n || throw(ArgumentError("joint histogram does not sum to samples"))
        for arm in 1:2
            l1,l2,f=state["logical_1_counts"][arm],state["logical_2_counts"][arm],state["failures"][arm]
            max(l1,l2) <= f <= min(n,l1+l2) || throw(ArgumentError("logical union count invariant failed"))
            marginal=sum(state["joint_counts"][mask+1] for mask in 0:3 if mask & (1 << (arm-1)) != 0)
            marginal == f || throw(ArgumentError("histogram marginal differs from arm failure count"))
        end
        isfinite(state["elapsed_seconds"]) && state["elapsed_seconds"] >= 0 || throw(ArgumentError("invalid elapsed time"))
    end
    return checkpoint
end

function csv_text(checkpoint)
    config=checkpoint["config"]; io=IOBuffer()
    columns=["channel","steps","Ly","Lz","p","shots","target_shots","seed","point_seed","complete","elapsed_seconds"]
    for arm in ARMS
        append!(columns,[arm*"_logical_1_failures",arm*"_logical_2_failures",arm*"_failures"])
    end
    append!(columns,["joint_"*lpad(mask,2,'0') for mask in 0:3])
    println(io,join(columns,','))
    for (index,point) in enumerate(point_grid(config))
        state=checkpoint["states"][index]
        row=Any[point.channel,point.L,point.L,point.L,point.p,state["completed_shots"],config["shots"],config["seed"],
            checkpoint["point_seeds_hex"][index],state["completed_shots"] == config["shots"],state["elapsed_seconds"]]
        for arm in 1:2
            append!(row,[state["logical_1_counts"][arm],state["logical_2_counts"][arm],state["failures"][arm]])
        end
        append!(row,state["joint_counts"])
        println(io,join(row,','))
    end
    return String(take!(io))
end

snapshot_name(checkpoint)=joinpath("checkpoints","checkpoint-"*lpad(checkpoint["generation"],6,'0')*"-"*checkpoint["payload_sha256"][1:12]*".toml")

function save_mirrors(directory,checkpoint)
    csv=csv_text(checkpoint)
    atomic_write(joinpath(directory,"raw.csv"),csv)
    states=checkpoint["states"]; target=checkpoint["config"]["shots"]
    metadata=Dict{String,Any}("format_version"=>1,"config"=>checkpoint["config"],"runtime"=>checkpoint["runtime"],
        "source_hashes"=>checkpoint["source_hashes"],"protocol"=>checkpoint["protocol"],
        "point_seeds_hex"=>checkpoint["point_seeds_hex"],"expected_points"=>length(states),
        "decoder_build_seconds_by_size"=>checkpoint["decoder_build_seconds_by_size"],
        "decoder_build_count_by_size"=>checkpoint["decoder_build_count_by_size"],
        "completed_points"=>count(s -> s["completed_shots"] == target,states),
        "completed_shots"=>sum(s["completed_shots"] for s in states),
        "complete"=>all(s["completed_shots"] == target for s in states),
        "raw_sha256"=>digest(csv),"checkpoint_sha256"=>digest(read(joinpath(directory,"checkpoint.toml"))),
        "current_snapshot"=>snapshot_name(checkpoint),"started_at_unix"=>checkpoint["started_at_unix"],
        "updated_at_unix"=>checkpoint["updated_at_unix"],"generation"=>checkpoint["generation"])
    atomic_write(joinpath(directory,"metadata.toml"),toml_text(metadata))
    return metadata
end

function commit_checkpoint!(directory,checkpoint)
    validate_states(checkpoint)
    pop!(checkpoint,"payload_sha256",nothing)
    checkpoint["generation"] += 1; checkpoint["updated_at_unix"]=time()
    checkpoint["payload_sha256"]=digest(toml_text(checkpoint))
    text=toml_text(checkpoint)
    snapshot=joinpath(directory,snapshot_name(checkpoint))
    mkpath(dirname(snapshot))
    if isfile(snapshot)
        read(snapshot,String) == text || error("immutable checkpoint snapshot collision")
    else
        atomic_write(snapshot,text)
    end
    # checkpoint.toml is authoritative; raw.csv and metadata.toml are mirrors
    # reconstructed on resume if a process stops between these atomic writes.
    atomic_write(joinpath(directory,"checkpoint.toml"),text)
    return save_mirrors(directory,checkpoint)
end

function load_checkpoint(directory,config,runtime,hashes)
    path=joinpath(directory,"checkpoint.toml")
    isfile(path) || throw(ArgumentError("resume checkpoint does not exist"))
    checkpoint=TOML.parsefile(path)
    claimed=pop!(checkpoint,"payload_sha256",nothing)
    claimed == digest(toml_text(checkpoint)) || throw(ArgumentError("checkpoint payload checksum mismatch"))
    checkpoint["payload_sha256"]=claimed
    for (key,expected) in (("format_version",1),("config",config),("runtime",runtime),
            ("source_hashes",hashes),("protocol",protocol_metadata(config["channel"])))
        get(checkpoint,key,nothing) == expected || throw(ArgumentError("resume mismatch: $key"))
    end
    validate_states(checkpoint)
    snapshot=joinpath(directory,snapshot_name(checkpoint))
    isfile(snapshot) && read(snapshot) == read(path) || throw(ArgumentError("immutable checkpoint snapshot mismatch"))
    return checkpoint
end

function dump_failure(directory,point,seed,shot,faults,error)
    record=Dict{String,Any}("channel"=>point.channel,"point_seed_hex"=>hexseed(seed),"shot_index"=>shot,
        "dimensions"=>[point.L,point.L,point.L],"p"=>point.p,"exception"=>sprint(showerror,error),
        "created_at_unix"=>time(),"fault_encoding"=>"1-based column-major linear indices of true bits; all omitted entries false")
    if faults !== nothing
        transfer=point.channel == "x_only" ? faults : faults.transfer
        if point.channel == "z_only"
            record["seed_length"]=length(faults.seed)
            record["seed_true_indices"]=findall(faults.seed)
        end
        record["transfer_shape"]=collect(size(transfer.z))
        record["transfer_x_true_indices"]=findall(vec(transfer.x))
        record["transfer_z_true_indices"]=findall(vec(transfer.z))
    end
    atomic_write(joinpath(directory,"failed-record.toml"),toml_text(record))
end

function two_outcomes(encoder,faults,decoder,channel)
    original=deepcopy(faults); results=Any[]
    for feedback in (false,true)
        result=channel == "x_only" ?
            replay_toric3d_h6_full_prefix_x(encoder,faults,decoder;feedback,record_history=false) :
            replay_toric3d_h6_full_prefix_z(encoder,faults,decoder;feedback,record_history=false)
        faults == original || error("full-prefix replay mutated the shared fault record")
        all(getproperty(result,k) isa Bool for k in (:logical_1,:logical_2,:failure)) || error("non-Boolean logical outcome")
        result.failure == (result.logical_1 | result.logical_2) || error("union flag differs from logical flags")
        push!(results,result)
    end
    return results
end

function run_study(options;io::IO=stdout)
    config=options.config; directory=options.output
    runtime=runtime_metadata(); hashes=source_hashes()
    if options.resume
        checkpoint=load_checkpoint(directory,config,runtime,hashes)
        metadata=save_mirrors(directory,checkpoint)
    else
        isfile(joinpath(directory,"checkpoint.toml")) && throw(ArgumentError("output already has a checkpoint; use --resume"))
        isdir(directory) && any(name -> name != "run.log",readdir(directory)) &&
            throw(ArgumentError("new output directory must be empty except for run.log"))
        mkpath(directory)
        checkpoint=fresh_checkpoint(config,runtime,hashes)
        metadata=commit_checkpoint!(directory,checkpoint)
    end
    metadata["complete"] && return metadata
    batches=0; cached_size=0; cached_p=NaN; context=nothing
    for (index,point) in enumerate(point_grid(config))
        state=checkpoint["states"][index]
        state["completed_shots"] == config["shots"] && continue
        seed=parse(UInt64,checkpoint["point_seeds_hex"][index];base=16)
        if point.L != cached_size || (point.channel == "x_only" && point.p != cached_p)
            build_seconds=@elapsed begin
                if point.L != cached_size
                    code=Toric3DEncodedSlice(point.L,point.L,point.L)
                    encoder=toric3d_encoded_encoder(code)
                else
                    code,encoder=context.code,context.encoder
                end
                decoder=point.channel == "x_only" ?
                    build_toric3d_x_feedback_decoders(code;p=point.p,settings=Toric3DDecoderSettings()) :
                    build_toric3d_z_feedback_decoders(code)
                context=(;code,encoder,decoder)
            end
            cached_size=point.L; cached_p=point.p
            checkpoint["decoder_build_seconds_by_size"][string(point.L)] += build_seconds
            checkpoint["decoder_build_count_by_size"][string(point.L)] += 1
            metadata=commit_checkpoint!(directory,checkpoint)
            println(io,"Built ",point.channel," L=",point.L," decoder context in ",round(build_seconds;digits=3)," s")
            flush(io)
        end
        encoder,decoder=context.encoder,context.decoder
        rng=resume_rng(seed,state["completed_shots"],encoder,point.p,point.channel)
        while state["completed_shots"] < config["shots"]
            count_shots=min(config["batch_size"],config["shots"]-state["completed_shots"])
            elapsed=@elapsed for _ in 1:count_shots
                shot=state["completed_shots"]+1; faults=nothing
                try
                    faults=sample_faults(rng,encoder,point.p,point.channel)
                    outcomes=two_outcomes(encoder,faults,decoder,point.channel)
                    for arm in 1:2
                        state["logical_1_counts"][arm] += outcomes[arm].logical_1
                        state["logical_2_counts"][arm] += outcomes[arm].logical_2
                        state["failures"][arm] += outcomes[arm].failure
                    end
                    state["joint_counts"][failure_mask([outcome.failure for outcome in outcomes])+1] += 1
                    state["completed_shots"] += 1
                catch error
                    error isa InterruptException && rethrow()
                    dump_failure(directory,point,seed,shot,faults,error)
                    rethrow()
                end
            end
            state["elapsed_seconds"] += elapsed
            metadata=commit_checkpoint!(directory,checkpoint)
            println(io,"H+6 ",point.channel," L=",point.L," p=",point.p," seed=",hexseed(seed),": ",state["completed_shots"],"/",config["shots"],
                " records; arm failures=",state["failures"],"; ",round(state["elapsed_seconds"];digits=3)," s")
            flush(io)
            batches+=1
            batches >= options.max_batches && return metadata
        end
    end
    return metadata
end

function main(arguments=ARGS;io::IO=stdout,error_io::IO=stderr)
    "--help" in arguments && (print(io,USAGE); return 0)
    try
        metadata=run_study(parse_arguments(arguments);io)
        println(io,"Saved ",metadata["completed_points"],"/",metadata["expected_points"]," completed points; ",metadata["completed_shots"]," paired records")
        return 0
    catch error
        error isa InterruptException && rethrow()
        println(error_io,"error: ",sprint(showerror,error)); print(error_io,USAGE)
        return 1
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(H6FullPrefixStudy.main())
end
