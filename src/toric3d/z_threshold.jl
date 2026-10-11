"""Exact final Z supports of independent seed-H and transfer-gate output faults."""
struct Toric3DZFaultModel
    code::Toric3DEncodedSlice
    encoder::Toric3DSliceEncoder
    matrix::SparseArrays.SparseMatrixCSC{UInt8,Int}
    locations::Vector{NamedTuple{(:kind,:step,:qubit),Tuple{Symbol,Int,Int}}}
end

"""
    build_toric3d_z_fault_model(code)

Include one noisy H on each seed edge before ideal compatibility measurements,
then every active output of the fixed transfer schedule. Every column is an
actual propagated data error. Sampling these columns preserves circuit-induced
correlations exactly; it does not substitute independent data-qubit noise.
"""
function build_toric3d_z_fault_model(code::Toric3DEncodedSlice)
    encoder = toric3d_encoded_encoder(code)
    n = data_qubit_count(code)
    seed = toric3d_slice_qubits(code.lattice,0)
    locations = [(kind=:seed,step=0,qubit=q) for q in seed]
    offsets = Int[]
    seen = falses(n); seen[seed] .= true
    for (step,operation) in enumerate(encoder.operations)
        operation.gate === :H && seen[only(operation.qubits)] &&
            error("Z-only suffix shortcut requires a fresh H target")
        push!(offsets,length(locations))
        for q in operation.qubits
            push!(locations,(kind=:gate,step=step,qubit=q))
            seen[q] = true
        end
    end
    supports = Vector{Vector{Int}}(undef,length(locations))
    suffix = [[q] for q in 1:n]
    for step in reverse(eachindex(encoder.operations))
        operation = encoder.operations[step]
        for (slot,q) in enumerate(operation.qubits)
            supports[offsets[step]+slot] = copy(suffix[q])
        end
        if operation.gate === :CNOT
            control,target = operation.qubits
            suffix[target] = sort!(symdiff(suffix[target],suffix[control]))
        elseif operation.gate === :H
            # This fresh qubit has no earlier active fault location. Its
            # pre-H Z image is never used, so no X suffix needs to be stored.
            nothing
        else
            error("unsupported gate in Z fault model")
        end
    end
    for (column,q) in enumerate(seed)
        suffix[q] == [q] || error("seed Z errors must remain on old controls")
        supports[column] = [q]
    end
    all(s -> 1 <= length(s) <= 3,supports) || error("Z fault support exceeded three qubits")
    rows,columns = Int[],Int[]
    for (column,support) in enumerate(supports), q in support
        push!(rows,q); push!(columns,column)
    end
    matrix = SparseArrays.sparse(rows,columns,ones(UInt8,length(rows)),n,length(supports))
    return Toric3DZFaultModel(code,encoder,matrix,locations)
end

"""Sample independent gate-fault bits and exactly propagate them in a batch."""
function sample_toric3d_z_frames(rng::Random.AbstractRNG,model::Toric3DZFaultModel;
        p::Real,shots::Integer)
    _validate_fault_rate(p,"p")
    count_shots = _positive_machine_int(shots,"shots")
    faults = Matrix{UInt8}(undef,size(model.matrix,2),count_shots)
    for index in eachindex(faults)
        faults[index] = rand(rng) < p
    end
    return BitMatrix(isodd.(model.matrix*faults))
end

struct Toric3DZThresholdDecoder
    code::Toric3DEncodedSlice
    matrix::SparseArrays.SparseMatrixCSC{UInt8,Int}
    logicals::SparseArrays.SparseMatrixCSC{UInt8,Int}
    matcher::PythonCall.Py
    version::String
end

function _toric3d_scipy_sparse(matrix::SparseArrays.SparseMatrixCSC)
    numpy = PythonCall.pyimport("numpy")
    scipy = PythonCall.pyimport("scipy.sparse")
    rows,columns,values = SparseArrays.findnz(matrix)
    return scipy.csc_matrix((numpy.asarray(UInt8.(values);dtype=numpy.uint8),
        (numpy.asarray(rows .- 1),numpy.asarray(columns .- 1)));shape=size(matrix))
end

"""Final unit-weight data-edge MWPM, with full data-qubit recovery vectors."""
function build_toric3d_z_threshold_decoder(code::Toric3DEncodedSlice)
    H = UInt8.(star_check_matrix(code;sparse=true))
    all(==(2),vec(sum(H;dims=1))) || error("encoded stars must form a closed graph")
    logicals = UInt8.(_check_matrix(toric3d_encoded_logical_x(code),data_qubit_count(code);sparse=true))
    pymatching = _load_pymatching()
    matcher = pymatching.Matching.from_check_matrix(_toric3d_scipy_sparse(H);
        weights=1.0,merge_strategy="smallest-weight",min_num_fault_ids=data_qubit_count(code))
    return Toric3DZThresholdDecoder(code,H,logicals,matcher,
        _toric3d_module_version(pymatching,"pymatching"))
end

function _toric3d_validate_z_batch(decoder::Toric3DZThresholdDecoder,syndromes,recoveries)
    size(recoveries) == (size(decoder.matrix,2),size(syndromes,2)) ||
        throw(DimensionMismatch("incorrect data-recovery batch dimensions"))
    BitMatrix(isodd.(decoder.matrix*Matrix{UInt8}(recoveries))) == syndromes ||
        error("invalid Z recovery: its star syndrome differs from the measurement")
    return recoveries
end

"""Decode only measured star syndromes; reject every invalid returned recovery."""
function decode_toric3d_z_batch(decoder::Toric3DZThresholdDecoder,
        syndromes::AbstractMatrix{Bool})
    size(syndromes,1) == size(decoder.matrix,1) ||
        throw(DimensionMismatch("incorrect star-syndrome batch dimensions"))
    all(iseven,vec(sum(syndromes;dims=1))) ||
        throw(ArgumentError("closed-graph star syndromes must have even parity"))
    n,shots = size(decoder.matrix,2),size(syndromes,2)
    shots == 0 && return falses(n,0)
    numpy = PythonCall.pyimport("numpy")
    input = numpy.asarray(Matrix{UInt8}(permutedims(syndromes));dtype=numpy.uint8)
    predicted = decoder.matcher.decode_batch(input)
    bits = _toric3d_py_bits(numpy,predicted,n*shots)
    recoveries = BitMatrix(reshape(bits,n,shots))
    return _toric3d_validate_z_batch(decoder,syndromes,recoveries)
end

"""Score either of two residual logical Z windings after validated recovery."""
function toric3d_z_threshold_counts(decoder::Toric3DZThresholdDecoder,
        frames::AbstractMatrix{Bool})
    size(frames,1) == data_qubit_count(decoder.code) ||
        throw(DimensionMismatch("incorrect data-frame dimensions"))
    syndrome = BitMatrix(isodd.(decoder.matrix*Matrix{UInt8}(frames)))
    recovery = decode_toric3d_z_batch(decoder,syndrome)
    residual = frames .⊻ recovery
    logical = isodd.(decoder.logicals*Matrix{UInt8}(residual))
    return (failures=count(vec(any(logical;dims=1))),
        logical_1_failures=count(@view logical[1,:]),
        logical_2_failures=count(@view logical[2,:]))
end

struct Toric3DZThresholdPoint
    dimensions::NTuple{3,Int}
    p::Float64
    shots::Int
    failures::Int
    logical_1_failures::Int
    logical_2_failures::Int
    elapsed_seconds::Float64
end

struct Toric3DZThresholdScan
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    shots::Int
    seed::Union{Nothing,Int}
    batch_size::Int
    points::Vector{Toric3DZThresholdPoint}
    rng_type::String
    decoder_version::String
    started_at_unix::Float64
end

function save_toric3d_z_threshold(scan::Toric3DZThresholdScan,directory::AbstractString)
    isempty(strip(directory)) && throw(ArgumentError("output directory must not be empty"))
    directory = abspath(directory); mkpath(directory)
    csv = joinpath(directory,"toric3d-z-threshold-raw.csv")
    _toric3d_atomic_write(csv) do io
        println(io,"steps,Ly,Lz,p,shots,failures,logical_1_failures,logical_2_failures,seed,decoder,elapsed_seconds")
        for point in scan.points
            entries = (point.dimensions...,point.p,point.shots,point.failures,
                point.logical_1_failures,point.logical_2_failures,something(scan.seed,"unspecified"),
                "unit_data_MWPM",point.elapsed_seconds)
            println(io,join(_toric3d_csv.(entries),','))
        end
    end
    metadata = joinpath(directory,"toric3d-z-threshold-metadata.toml")
    settings = Dict{String,Any}(
        "format_version"=>1,"dimensions_steps_Ly_Lz"=>[collect(d) for d in scan.dimensions],
        "error_rates"=>scan.error_rates,"shots_per_point"=>scan.shots,
        "seed"=>something(scan.seed,"unspecified"),"rng_type"=>scan.rng_type,
        "batch_size"=>scan.batch_size,"julia_version"=>string(VERSION),
        "started_at_unix"=>scan.started_at_unix,"completed_points"=>length(scan.points),
        "expected_points"=>length(scan.dimensions)*length(scan.error_rates),
        "complete"=>length(scan.points)==length(scan.dimensions)*length(scan.error_rates),
        "logical_qubits"=>2,"logical_z_distance"=>"2 min(Ly,Lz)",
        "qubit_count"=>"3 (steps+1) Ly Lz",
        "seed_preparation"=>"noisy H on each of 3 Ly Lz initial qubits, then ideal Q measurements and ideal X recovery; logical |++>",
        "transfer_circuit"=>"same diagonal seven-gate triad schedule as slice experiment",
        "fault_locations"=>"3 Ly Lz seed H outputs + 13 steps Ly Lz transfer gate outputs",
        "noise"=>"independent Bernoulli-p Z faults on every participating H/CNOT output",
        "ideal_operations"=>["initial product zeros","measurements","corrections","idles"],
        "sampling"=>"exact sparse GF(2) propagation of independent gate-fault bits; no iid data-noise replacement",
        "decoder"=>"unit-weight data-edge PyMatching MWPM; full recovery checked against measured syndrome",
        "decoder_version"=>scan.decoder_version,
        "decoder_limitations"=>"data-edge weights ignore known local circuit-fault correlations; decoder-specific threshold",
        "feedback"=>"pure-Z faults do not trigger Q measurements; active-slice feedback is identically inactive",
        "failure"=>"either residual transverse logical-Z winding after final recovery; equals |++> state failure in this Z-only study",
        "scope"=>"encoded diagonal slab with noisy measurement-prepared seed; not a full X/Z threshold or the earlier k=0 defect-density experiment")
    _toric3d_atomic_write(metadata) do io
        TOML.print(io,settings;sorted=true)
    end
    return (csv=csv,metadata=metadata)
end

"""
    scan_toric3d_z_threshold(rng; sizes=[5,7,9], error_rates=..., shots=1000,
        seed=nothing, batch_size=256, output_dir=nothing, progress_io=nothing)

Simulate logical-|++> preparation with noisy seed H gates, ideal seed
measurements, unchanged transfer-gate Z noise, and final unit-data-weight MWPM.
The supplied RNG controls all independent fault draws; seed is provenance.
Save completed points incrementally when output_dir is given.
"""
function scan_toric3d_z_threshold(rng::Random.AbstractRNG;
        sizes=[5,7,9],error_rates=[.002,.004,.006,.008,.010,.015,.020],
        shots::Integer=1000,seed=nothing,batch_size::Integer=256,
        output_dir=nothing,progress_io=nothing)
    dimensions,rates = _toric3d_sizes(sizes),_scan_error_rates(error_rates)
    count_shots = _positive_machine_int(shots,"shots")
    batch = _positive_machine_int(batch_size,"batch_size")
    progress_io === nothing || progress_io isa IO || throw(ArgumentError("progress_io must be an IO"))
    output_dir === nothing || output_dir isa AbstractString || throw(ArgumentError("output_dir must be a path"))
    scan = Toric3DZThresholdScan(dimensions,rates,count_shots,_seed_metadata(seed),batch,
        Toric3DZThresholdPoint[],string(typeof(rng)),_toric3d_module_version(_load_pymatching(),"pymatching"),time())
    for dims in dimensions
        code = Toric3DEncodedSlice(dims...)
        model,decoder = build_toric3d_z_fault_model(code),build_toric3d_z_threshold_decoder(code)
        for p in rates
            started = time(); failures = logical1 = logical2 = 0
            for first in 1:batch:count_shots
                size_batch = min(batch,count_shots-first+1)
                frames = sample_toric3d_z_frames(rng,model;p,shots=size_batch)
                try
                    result = toric3d_z_threshold_counts(decoder,frames)
                    failures += result.failures
                    logical1 += result.logical_1_failures
                    logical2 += result.logical_2_failures
                catch error
                    throw(ErrorException("Z threshold recovery failed at dimensions=$dims,p=$p,seed=$(scan.seed),batch_start=$first: "*sprint(showerror,error)))
                end
            end
            point = Toric3DZThresholdPoint(dims,p,count_shots,failures,logical1,logical2,time()-started)
            push!(scan.points,point)
            output_dir === nothing || save_toric3d_z_threshold(scan,output_dir)
            if progress_io !== nothing
                println(progress_io,"Z threshold ",dims," p=",p,": ",failures,"/",count_shots,
                    " = ",round(failures/count_shots;digits=5)," in ",round(point.elapsed_seconds;digits=2)," s")
                flush(progress_io)
            end
        end
    end
    return scan
end
