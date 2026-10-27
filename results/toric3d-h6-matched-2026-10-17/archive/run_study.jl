module H6MatchedStudy

using TopoNoise, Random, LinearAlgebra
import TOML, SHA, PythonCall

const ROOT=abspath(joinpath(@__DIR__,"runtime"))
const PYTHON="/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python"
const ARMS=("final_only","one_step","two_slice")
const USAGE="""
Usage: julia --project=. --startup-file=no run_study.jl --channel x_only|z_only --output DIR --sizes 3,7,11 --rates 0,.005,.01 --shots N --seed N [--batch-size 25] [--resume] [--max-batches N]

Matched legacy H+6 final-only, one-step, and complete two-slice recovery, T=L.
Each channel uses one immutable fault record per triple and a shared final decoder.
X starts from ideal |00>. Z samples noisy seed-H preparation and applies the
same complete r=0 seed syndrome repair to all three arms, retaining logical errors.
One-step timing follows the existing channel API; two-slice acts after 1:T.
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

failure_mask(flags)=sum(Int(flags[i]) << (i-1) for i in 1:3)
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

function source_hashes(root=ROOT,runner=@__FILE__)
    isdir(joinpath(root,"src")) || throw(ArgumentError("frozen runtime src directory missing"))
    paths=["Project.toml","Manifest.toml"]
    for (directory,_,files) in walkdir(joinpath(root,"src")), file in files
        push!(paths,relpath(joinpath(directory,file),root))
    end
    hashes=Dict("runtime/"*path=>digest(read(joinpath(root,path))) for path in sort!(paths))
    hashes["run_study.jl"]=digest(read(runner))
    return hashes
end

function runtime_metadata()
    realpath(dirname(Base.active_project())) == realpath(ROOT) ||
        throw(ArgumentError("activate this archive/runtime project"))
    realpath(pathof(TopoNoise)) == realpath(joinpath(ROOT,"src","TopoNoise.jl")) ||
        throw(ArgumentError("TopoNoise must be loaded from the frozen runtime"))
    Threads.nthreads() == 1 && BLAS.get_num_threads() == 1 ||
        throw(ArgumentError("one Julia and BLAS thread required"))
    thread_keys=("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    all(get(ENV,key,"") == "1" for key in thread_keys) ||
        throw(ArgumentError("all four thread environment variables must equal 1"))
    get(ENV,"JULIA_CONDAPKG_BACKEND","") == "Null" || throw(ArgumentError("JULIA_CONDAPKG_BACKEND=Null required"))
    executable=PythonCall.pyconvert(String,PythonCall.pyimport("sys").executable)
    requested=get(ENV,"JULIA_PYTHONCALL_EXE","")
    !isempty(requested) && realpath(executable) == realpath(requested) == realpath(PYTHON) ||
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
        "initial_state"=>channel == "x_only" ? "ideal encoded logical |00>" : "noisy seed-H preparation of logical |++> followed by common r=0 repair",
        "seed_preparation"=>channel == "x_only" ? "ideal; zero seed fault locations" : "3 L^2 noisy seed-H outputs; identical complete r=0 syndrome repair supplied to every arm; residual logical errors retained",
        "common_z_seed_repair"=>"syndrome=isodd.(family.seed.matrix*UInt8.(raw.seed)); recovery=decode_toric3d_z_feedback(family.seed,syndrome); corrected.seed=raw.seed XOR recovery; transfer faults unchanged; corrected seed syndrome zero",
        "final_only"=>"no intermediate correction after common preparation; exactly one shared final recovery and logical outcome",
        "one_step"=>channel == "x_only" ? "existing replay_toric3d_x_feedback with routing=:one_step and two.final" : "existing replay_toric3d_h6_z_one_step on the common corrected record; redundant internal seed correction is zero; retiring-sheet endpoint corrections after transfers 1:T-1",
        "two_slice"=>"complete physical recovery on newest two data slices after EVERY transfer 1:T, including T; complete-prefix syndrome zero; exactly one shared final recovery and logical outcome",
        "feedback_timing"=>"common Z r=0 repair; one-step uses existing channel timing; complete two-slice after transfers 1:T inclusive",
        "fault_locations"=>channel == "x_only" ? "13 L^3 participating transfer outputs" : "3 L^2 seed-H outputs plus 13 L^3 participating transfer outputs",
        "noise"=>"independent Bernoulli-p faults in the selected channel at original recorded circuit locations",
        "sampler"=>channel == "x_only" ? "sample_toric3d_faults(p_x=p,p_z=0)" : "sample_toric3d_diagonal_z_faults: original noisy seed first, then participating transfer outputs in schedule order",
        "shared_final_decoder"=>channel == "x_only" ? "two.final; full-plaquette BP+OSD with declared p prior; syndrome-validated recovery" : "family.prefixes[end]; unit-weight data-edge PyMatching MWPM; syndrome-validated recovery",
        "decoder_settings"=>channel == "x_only" ? Dict(string(k)=>v for (k,v) in pairs(settings)) : Dict("weights"=>1.0,"merge_strategy"=>"smallest-weight"),
        "z_weights"=>"unit data-edge weights","failure"=>"either residual logical parity after final recovery",
        "joint_histogram_encoding"=>"mask = final_only + 2*one_step + 4*two_slice; each term is its state-failure Boolean; joint_00 through joint_07 are decimal masks",
        "point_order"=>"sizes outer loop, rates inner loop; both strictly increasing",
        "point_seed_generation"=>"MersenneTwister(integer base seed); repeated scalar rand(rng,UInt64); discard duplicate values; 16 lowercase hexadecimal digits",
        "point_rng"=>"Random.Xoshiro(point_seed)",
        "resume"=>"recreate point Xoshiro and sample/discard completed ORIGINAL records without decoding; reject corrupt authoritative checkpoint; rebuild raw CSV and metadata mirrors")
end

function fresh_checkpoint(config,runtime,hashes)
    count=length(point_grid(config))
    states=[Dict{String,Any}("completed_shots"=>0,"logical_1_counts"=>zeros(Int,3),
        "logical_2_counts"=>zeros(Int,3),"failures"=>zeros(Int,3),"joint_counts"=>zeros(Int,8),
        "elapsed_seconds"=>0.0) for _ in 1:count]
    return Dict{String,Any}("format_version"=>1,"config"=>config,"runtime"=>runtime,
        "source_hashes"=>hashes,"protocol"=>protocol_metadata(config["channel"]),"point_seeds_hex"=>hexseed.(point_seeds(config["seed"],count)),
        "states"=>states,"generation"=>0,"started_at_unix"=>time(),"updated_at_unix"=>time(),
        "decoder_build_seconds_by_size"=>Dict(string(L)=>0.0 for L in config["sizes"]),
        "decoder_build_count_by_size"=>Dict(string(L)=>0 for L in config["sizes"]))
end

function validate_states(checkpoint)
    config=checkpoint["config"]; states=checkpoint["states"]
    generation=checkpoint["generation"]
    generation isa Integer && !(generation isa Bool) && generation >= 0 ||
        throw(ArgumentError("invalid checkpoint generation"))
    started=checkpoint["started_at_unix"]; updated=checkpoint["updated_at_unix"]
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x >= 0,(started,updated)) && updated >= started ||
        throw(ArgumentError("invalid checkpoint timestamps"))
    length(states) == length(point_grid(config)) || throw(ArgumentError("checkpoint grid length mismatch"))
    checkpoint["point_seeds_hex"] == hexseed.(point_seeds(config["seed"],length(states))) ||
        throw(ArgumentError("checkpoint point seeds differ from exact base-seed regeneration"))
    for L in config["sizes"]
        seconds=checkpoint["decoder_build_seconds_by_size"][string(L)]
        builds=checkpoint["decoder_build_count_by_size"][string(L)]
        seconds isa Real && !(seconds isa Bool) && isfinite(seconds) && seconds >= 0 &&
            builds isa Integer && !(builds isa Bool) && builds >= 0 ||
            throw(ArgumentError("invalid decoder build timing"))
    end
    unfinished=false
    for state in states
        n=state["completed_shots"]
        n isa Integer && !(n isa Bool) && 0 <= n <= config["shots"] || throw(ArgumentError("invalid completed sample count"))
        unfinished && n != 0 && throw(ArgumentError("checkpoint points violate sequential sampling order"))
        n < config["shots"] && (unfinished=true)
        for (key,size) in (("logical_1_counts",3),("logical_2_counts",3),("failures",3),("joint_counts",8))
            counts=state[key]
            length(counts) == size && all(x -> x isa Integer && !(x isa Bool) && 0 <= x <= n,counts) ||
                throw(ArgumentError("invalid $key"))
        end
        sum(state["joint_counts"]) == n || throw(ArgumentError("joint histogram does not sum to samples"))
        for arm in 1:3
            l1,l2,f=state["logical_1_counts"][arm],state["logical_2_counts"][arm],state["failures"][arm]
            max(l1,l2) <= f <= min(n,l1+l2) || throw(ArgumentError("logical union count invariant failed"))
            marginal=sum(state["joint_counts"][mask+1] for mask in 0:7 if mask & (1 << (arm-1)) != 0)
            marginal == f || throw(ArgumentError("histogram marginal differs from arm failure count"))
        end
        seconds=state["elapsed_seconds"]
        seconds isa Real && !(seconds isa Bool) && isfinite(seconds) && seconds >= 0 ||
            throw(ArgumentError("invalid elapsed time"))
    end
    return checkpoint
end

function csv_text(checkpoint)
    config=checkpoint["config"]; io=IOBuffer()
    columns=["channel","steps","Ly","Lz","p","shots","target_shots","seed","point_seed","complete","elapsed_seconds"]
    for arm in ARMS
        append!(columns,[arm*"_logical_1_failures",arm*"_logical_2_failures",arm*"_failures"])
    end
    append!(columns,["joint_"*lpad(mask,2,'0') for mask in 0:7])
    println(io,join(columns,','))
    for (index,point) in enumerate(point_grid(config))
        state=checkpoint["states"][index]
        row=Any[point.channel,point.L,point.L,point.L,point.p,state["completed_shots"],config["shots"],config["seed"],
            checkpoint["point_seeds_hex"][index],state["completed_shots"] == config["shots"],state["elapsed_seconds"]]
        for arm in 1:3
            append!(row,[state["logical_1_counts"][arm],state["logical_2_counts"][arm],state["failures"][arm]])
        end
        append!(row,state["joint_counts"])
        println(io,join(row,','))
    end
    return String(take!(io))
end


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
        "previous_checkpoint"=>"previous-checkpoint.toml","started_at_unix"=>checkpoint["started_at_unix"],
        "updated_at_unix"=>checkpoint["updated_at_unix"],"generation"=>checkpoint["generation"])
    atomic_write(joinpath(directory,"metadata.toml"),toml_text(metadata))
    return metadata
end

function read_checked_checkpoint(path)
    isfile(path) || throw(ArgumentError("checkpoint does not exist: $path"))
    checkpoint=try
        TOML.parsefile(path)
    catch error
        throw(ArgumentError("checkpoint is not valid TOML: "*sprint(showerror,error)))
    end
    claimed=pop!(checkpoint,"payload_sha256",nothing)
    claimed isa String && claimed == digest(toml_text(checkpoint)) ||
        throw(ArgumentError("checkpoint payload checksum mismatch: $path"))
    checkpoint["payload_sha256"]=claimed
    return checkpoint
end

function commit_checkpoint!(directory,checkpoint)
    validate_states(checkpoint)
    current=joinpath(directory,"checkpoint.toml")
    if isfile(current)
        saved=read_checked_checkpoint(current)
        get(checkpoint,"payload_sha256",nothing) == saved["payload_sha256"] ||
            throw(ArgumentError("authoritative checkpoint changed during this run"))
        # Backup is itself checksummed. If interrupted here the primary remains authoritative.
        atomic_write(joinpath(directory,"previous-checkpoint.toml"),read(current,String))
    end
    pop!(checkpoint,"payload_sha256",nothing)
    checkpoint["generation"] += 1; checkpoint["updated_at_unix"]=time()
    checkpoint["payload_sha256"]=digest(toml_text(checkpoint))
    atomic_write(current,toml_text(checkpoint))
    return save_mirrors(directory,checkpoint)
end

function load_checkpoint(directory,config,runtime,hashes)
    checkpoint=read_checked_checkpoint(joinpath(directory,"checkpoint.toml"))
    for (key,expected) in (("format_version",1),("config",config),("runtime",runtime),
            ("source_hashes",hashes),("protocol",protocol_metadata(config["channel"])))
        get(checkpoint,key,nothing) == expected || throw(ArgumentError("resume mismatch: $key"))
    end
    validate_states(checkpoint)
    backup=joinpath(directory,"previous-checkpoint.toml")
    if isfile(backup)
        previous=read_checked_checkpoint(backup)
        for key in ("format_version","config","runtime","source_hashes","protocol","point_seeds_hex")
            get(previous,key,nothing) == checkpoint[key] || throw(ArgumentError("previous checkpoint mismatch: $key"))
        end
        validate_states(previous)
        0 <= checkpoint["generation"]-previous["generation"] <= 1 ||
            throw(ArgumentError("previous checkpoint generation mismatch"))
    end
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

function build_context(L,p,channel; previous=nothing)
    channel in ("x_only","z_only") || throw(ArgumentError("unsupported noise channel"))
    reusable=previous !== nothing && previous.code.lattice.steps == L && previous.channel == channel
    code=reusable ? previous.code : Toric3DEncodedSlice(L,L,L)
    encoder=reusable ? previous.encoder : toric3d_encoded_encoder(code)
    if channel == "x_only"
        slice=reusable ? previous.slice : build_toric3d_slice_decoder(code.lattice)
        two=build_toric3d_h6_two_slice_x_decoder(code;p,settings=Toric3DDecoderSettings())
        return (;code,encoder,slice,two,channel)
    end
    family=reusable ? previous.family : build_toric3d_z_feedback_decoders(code)
    one=reusable ? previous.one : build_toric3d_h6_z_one_step_decoder(family)
    return (;code,encoder,family,one,channel)
end

"""Repair only the noisy seed syndrome; never erase its residual logical sector."""
function common_seed_record(raw,family)
    family.seed.window.qubits == toric3d_slice_qubits(family.code.lattice,0) ||
        error("seed decoder qubit ordering differs from fault record")
    syndrome=BitVector(isodd.(family.seed.matrix*UInt8.(raw.seed)))
    recovery=decode_toric3d_z_feedback(family.seed,syndrome)
    corrected=Toric3DDiagonalZFaultRecord(raw.seed .⊻ recovery,raw.transfer)
    any(isodd.(family.seed.matrix*UInt8.(corrected.seed))) && error("common seed repair left a syndrome")
    return corrected,syndrome,recovery
end

function three_outcomes(context,raw,channel;record_history=false)
    context.channel == channel || throw(ArgumentError("decoder context channel mismatch"))
    original=deepcopy(raw)
    corrected=channel == "z_only" ? first(common_seed_record(raw,context.family)) : raw
    corrected_original=channel == "z_only" ? deepcopy(corrected) : original
    results=Any[]
    for arm in 1:3
        result=if channel == "x_only"
            arm < 3 ? replay_toric3d_x_feedback(context.encoder,raw,context.slice,context.two.final;
                feedback=arm == 2,routing=:one_step,record_history) :
                replay_toric3d_h6_two_slice_x(context.encoder,raw,context.two;feedback=true,record_history)
        else
            arm == 2 ? replay_toric3d_h6_z_one_step(context.encoder,corrected,context.one;record_history) :
                replay_toric3d_h6_two_slice_z(context.encoder,corrected,context.family;feedback=arm == 3,record_history)
        end
        raw == original || error("replay mutated the original shared fault record")
        corrected == corrected_original || error("replay mutated the common corrected fault record")
        all(getproperty(result,k) isa Bool for k in (:logical_1,:logical_2,:failure)) || error("non-Boolean logical outcome")
        result.failure == (result.logical_1 | result.logical_2) || error("union flag differs from logical flags")
        push!(results,result)
    end
    return results
end

function accumulate!(state,outcomes)
    length(outcomes) == 3 || throw(ArgumentError("three outcomes required per shared fault record"))
    for result in outcomes
        all(getproperty(result,k) isa Bool for k in (:logical_1,:logical_2,:failure)) &&
            result.failure == (result.logical_1 | result.logical_2) || throw(ArgumentError("invalid logical outcome"))
    end
    for arm in 1:3
        state["logical_1_counts"][arm] += outcomes[arm].logical_1
        state["logical_2_counts"][arm] += outcomes[arm].logical_2
        state["failures"][arm] += outcomes[arm].failure
    end
    state["joint_counts"][failure_mask([outcome.failure for outcome in outcomes])+1] += 1
    state["completed_shots"] += 1
    return state
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
            build_seconds=@elapsed context=build_context(point.L,point.p,point.channel;previous=context)
            cached_size=point.L; cached_p=point.p
            checkpoint["decoder_build_seconds_by_size"][string(point.L)] += build_seconds
            checkpoint["decoder_build_count_by_size"][string(point.L)] += 1
            metadata=commit_checkpoint!(directory,checkpoint)
            println(io,"Built ",point.channel," L=",point.L," decoder context in ",round(build_seconds;digits=3)," s")
            flush(io)
        end
        encoder=context.encoder
        rng=resume_rng(seed,state["completed_shots"],encoder,point.p,point.channel)
        while state["completed_shots"] < config["shots"]
            count_shots=min(config["batch_size"],config["shots"]-state["completed_shots"])
            elapsed=@elapsed for _ in 1:count_shots
                shot=state["completed_shots"]+1; faults=nothing
                try
                    faults=sample_faults(rng,encoder,point.p,point.channel)
                    outcomes=three_outcomes(context,faults,point.channel)
                    accumulate!(state,outcomes)
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
        println(io,"Saved ",metadata["completed_points"],"/",metadata["expected_points"]," completed points; ",metadata["completed_shots"]," shared three-arm records")
        return 0
    catch error
        error isa InterruptException && rethrow()
        println(error_io,"error: ",sprint(showerror,error)); print(error_io,USAGE)
        return 1
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(H6MatchedStudy.main())
end
