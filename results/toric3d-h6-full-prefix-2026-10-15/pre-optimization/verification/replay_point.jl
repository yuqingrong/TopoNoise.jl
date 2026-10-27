"""Read-only point replay: independent RNG, propagation, geometry and scoring."""
module H6FullPrefixIndependentReplay

using TopoNoise, Random, LinearAlgebra, SparseArrays
import TOML, SHA, PythonCall

const ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
const ARMS = ("final_only", "full_prefix")
const USAGE = "Usage: replay_point.jl INPUT_DIR L P OUTPUT_TOML"

digest(bytes) = bytes2hex(SHA.sha256(bytes))
parity(matrix, bits) = BitVector(isodd.(matrix * Vector{UInt8}(bits)))
measured(checks, bits) = BitVector(isodd(count(bits[check])) for check in checks)

function check!(counts, category, condition, message)
    condition || error("independent replay [$category]: $message")
    counts[category] = get(counts, category, 0) + 1
    return nothing
end

function canonical_toml(data)
    io = IOBuffer()
    TOML.print(io, data; sorted=true)
    return String(take!(io))
end

function regenerate_seeds(base_seed, number)
    stream = MersenneTwister(base_seed)
    unique_seeds = Set{UInt64}()
    seeds = UInt64[]
    while length(seeds) < number
        seed = rand(stream, UInt64) # Scalar draws: array draws use a different stream.
        seed in unique_seeds && continue
        push!(unique_seeds, seed)
        push!(seeds, seed)
    end
    return seeds
end

"""Reproduce each documented scalar RNG draw without calling a production sampler."""
function regenerate_faults(rng, encoder, p, channel)
    lattice = encoder.code
    x, z = falses(2, length(encoder.operations)), falses(2, length(encoder.operations))
    seed = channel == "z_only" ? BitVector(rand(rng) < p for _ in 1:3lattice.Ly*lattice.Lz) : BitVector()
    for (operation_index, operation) in enumerate(encoder.operations), slot in eachindex(operation.qubits)
        if channel == "x_only"
            x[slot, operation_index] = rand(rng) < p
            # The established X-only sampler still consumes the p_z=0 draw.
            z[slot, operation_index] = rand(rng) < 0.0
        else
            z[slot, operation_index] = rand(rng) < p
        end
    end
    transfer = Toric3DFaultRecord(lattice.steps, lattice.Ly, lattice.Lz, x, z)
    return channel == "x_only" ? transfer : Toric3DDiagonalZFaultRecord(seed, transfer)
end

"""Apply the elementary Clifford rules directly, then inject output faults."""
function propagate_independently!(frame, encoder, transfer, first, last)
    for operation_index in first:last
        operation = encoder.operations[operation_index]
        if operation.gate === :H
            q = only(operation.qubits)
            frame.x[q], frame.z[q] = frame.z[q], frame.x[q]
        elseif operation.gate === :CNOT
            control, target = operation.qubits
            frame.x[target] ⊻= frame.x[control]
            frame.z[control] ⊻= frame.z[target]
        else
            error("independent replay: unexpected gate $(operation.gate)")
        end
        for (slot, q) in enumerate(operation.qubits)
            frame.x[q] ⊻= transfer.x[slot, operation_index]
            frame.z[q] ⊻= transfer.z[slot, operation_index]
        end
    end
    return frame
end

function check_runtime!(counts, saved)
    check!(counts, "runtime", string(VERSION) == saved["julia_version"], "Julia version differs")
    check!(counts, "runtime", Threads.nthreads() == BLAS.get_num_threads() == 1,
        "Julia and BLAS must each use one thread")
    check!(counts, "runtime", saved["julia_threads"] == saved["blas_threads"] == 1,
        "saved Julia/BLAS thread declaration differs")
    for name in ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        check!(counts, "runtime", get(ENV, name, "") == "1", "$name must equal 1")
    end
    check!(counts, "runtime", get(ENV, "JULIA_CONDAPKG_BACKEND", "") == "Null",
        "JULIA_CONDAPKG_BACKEND must equal Null")
    python = PythonCall.pyconvert(String, PythonCall.pyimport("sys").executable)
    check!(counts, "runtime", realpath(python) == realpath(saved["python_executable"]) ==
        realpath(ENV["JULIA_PYTHONCALL_EXE"]), "Python executable differs")
    versions = PythonCall.pyimport("importlib.metadata")
    check!(counts, "runtime", Set(keys(saved["packages"])) == Set(("numpy", "scipy", "ldpc", "pymatching")),
        "Python package manifest is incomplete or unexpected")
    python_version = PythonCall.pyconvert(String, PythonCall.pyimport("sys").version)
    check!(counts, "runtime", python_version == saved["python_version"], "Python version differs")
    environment_keys = ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "JULIA_CONDAPKG_BACKEND", "JULIA_PYTHONCALL_EXE")
    check!(counts, "runtime", Set(keys(saved["environment"])) == Set(environment_keys),
        "runtime environment manifest differs")
    for name in environment_keys
        check!(counts, "runtime", saved["environment"][name] == get(ENV, name, ""),
            "saved runtime environment differs for $name")
    end
    for (name, version) in saved["packages"]
        observed = PythonCall.pyconvert(String, versions.version(name))
        check!(counts, "runtime", observed == version, "$name version differs")
    end
    return Dict("julia_version"=>string(VERSION), "python_executable"=>python,
        "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads(),
        "python_version"=>python_version, "packages"=>saved["packages"],
        "environment"=>saved["environment"])
end

function incidence(checks, qubits)
    columns = Dict(q=>i for (i,q) in enumerate(qubits))
    rows, cols = Int[], Int[]
    for (i, check) in enumerate(checks), q in check
        column = get(columns, q, 0)
        column == 0 && continue
        push!(rows, i); push!(cols, column)
    end
    return sparse(rows, cols, ones(Int, length(rows)), length(checks), length(qubits))
end

function independent_geometry(code, channel)
    physical_checks = channel == "x_only" ? plaquette_checks : star_checks
    full_matrix = incidence(physical_checks(code), 1:data_qubit_count(code))
    logical_supports = channel == "x_only" ? toric3d_encoded_logical_z(code) : toric3d_encoded_logical_x(code)
    logical_matrix = incidence(logical_supports, 1:data_qubit_count(code))
    prefixes = map(1:code.lattice.steps) do r
        # Construct a fresh prefix and coordinate map without production prefix
        # builders, decoder matrices or history fields.
        prefix = Toric3DEncodedSlice(r, code.lattice.Ly, code.lattice.Lz)
        qubits = [data_qubit_index(code, data_qubit_coordinate(prefix, q)...)
            for q in 1:data_qubit_count(prefix)]
        checks = [qubits[check] for check in physical_checks(prefix)]
        future = setdiff(collect(1:data_qubit_count(code)), qubits)
        frontier = channel == "x_only" ?
            [qubits[check] for check in toric3d_slice_checks(prefix.lattice, r)] : Vector{Int}[]
        (;qubits, checks, matrix=incidence(checks, qubits), frontier, future)
    end
    return (;full_matrix, logical_matrix, prefixes)
end

function check_final!(counts, result, geometry, channel)
    bits = channel == "x_only" ? result.frame.x : result.frame.z
    other = channel == "x_only" ? result.frame.z : result.frame.x
    matrix = geometry.full_matrix
    check!(counts, "final", length(bits) == length(result.recovery) == length(result.residual) == size(matrix, 2),
        "final data dimensions differ")
    check!(counts, "final", !any(other), "single-channel replay generated the other Pauli channel")
    check!(counts, "final", parity(matrix, bits) == result.syndrome, "final physical syndrome differs")
    check!(counts, "final", parity(matrix, result.recovery) == result.syndrome, "final recovery syndrome differs")
    check!(counts, "final", result.residual == bits .⊻ result.recovery, "final residual is not frame XOR recovery")
    check!(counts, "final", !any(parity(matrix, result.residual)), "final residual is not syndrome-free")
    logical = parity(geometry.logical_matrix, result.residual)
    check!(counts, "final", (result.logical_1, result.logical_2) == Tuple(logical), "logical parities differ")
    check!(counts, "final", result.failure == any(logical), "logical union differs")
end

function check_history!(counts, result, encoder, faults, code, geometry, channel; feedback=true)
    steps = code.lattice.steps
    check!(counts, "history", length(result.history) == (feedback ? steps : 0),
        "history must contain exactly T prefixes with feedback, none without; no r=0 seed entry")
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    if channel == "z_only"
        seed_qubits = [q for q in 1:data_qubit_count(code) if data_qubit_coordinate(code, q)[2] == 0]
        frame.z[seed_qubits] .= faults.seed
    end
    transfer = channel == "x_only" ? faults : faults.transfer
    first_operation = 1
    for r in 1:steps
        propagate_independently!(frame, encoder, transfer, first_operation, encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        feedback || continue
        entry, prefix = result.history[r], geometry.prefixes[r]
        bits = channel == "x_only" ? frame.x : frame.z
        before = measured(prefix.checks, bits)
        reported = channel == "x_only" ? entry.pre_syndrome : entry.syndrome
        check!(counts, "history", entry.checkpoint == r, "prefix order differs or final transfer omitted")
        check!(counts, "history", entry.qubits == prefix.qubits, "recovery columns differ from complete Q0:Qr prefix")
        check!(counts, "history", length(entry.recovery) == length(prefix.qubits), "recovery length differs from prefix")
        check!(counts, "history", reported == before, "measurement differs from independently propagated frame")
        check!(counts, "history", parity(prefix.matrix, entry.recovery) == before,
            "recovery does not reproduce every physical prefix check")
        support = prefix.qubits[findall(entry.recovery)]
        check!(counts, "history", all(0 <= data_qubit_coordinate(code, q)[2] <= r for q in support),
            "recovery leaves completed Q0:Qr prefix")
        future_before = copy(bits[prefix.future])
        bits[prefix.qubits] .⊻= entry.recovery
        future_after = bits[prefix.future]
        check!(counts, "history", future_before == future_after, "recovery changed future data")
        post = measured(prefix.checks, bits)
        check!(counts, "history", !any(post), "complete-prefix syndrome remains nonzero")
        check!(counts, "history", entry.post_syndrome == post, "reported complete-prefix post-syndrome differs")
        if channel == "x_only"
            frontier = measured(prefix.frontier, bits)
            check!(counts, "history", !any(frontier) && entry.frontier_syndrome == frontier, "X frontier is not clear")
            check!(counts, "history", entry.future_qubits == prefix.future, "X future data indices differ")
            check!(counts, "history", entry.future_before == future_before && entry.future_after == future_after,
                "X future-frame history differs")
        else
            check!(counts, "history", entry.kind === :prefix, "Z history includes a seed or non-prefix correction")
            check!(counts, "history", entry.support == support, "Z physical recovery support differs")
        end
    end
    check!(counts, "history", frame.x == result.frame.x && frame.z == result.frame.z,
        "independent history propagation differs from pre-final frame")
    if feedback
        check!(counts, "history", !any(result.syndrome), "last complete prefix did not clear final syndrome")
        check!(counts, "history", !any(result.recovery), "final zero-syndrome recovery is nonzero")
    end
end

function check_protocol!(counts, protocol, channel)
    expected = Dict(
        "circuit"=>"legacy H+6-CNOT plaquette transfer", "geometry"=>"T=L; Ly=Lz=L",
        "channel"=>channel, "arm_order"=>collect(ARMS),
        "initial_state"=>channel == "x_only" ? "ideal encoded logical |00>" : "noisy seed-H preparation of logical |++>",
        "seed_preparation"=>channel == "x_only" ? "ideal; zero seed fault locations" : "3 L^2 noisy seed-H outputs; no separate r=0 seed correction in either arm",
        "final_only"=>"no intermediate physical correction; exactly one final recovery and logical outcome",
        "full_prefix"=>"complete physical recovery on all available data slices Q0..Qr after EVERY transfer 1:T, including the last; complete-prefix syndrome zero; exactly one final recovery and logical outcome afterward",
        "feedback_timing"=>"after transfers 1:T inclusive; no separate r=0 correction",
        "fault_locations"=>channel == "x_only" ? "13 L^3 participating transfer outputs" : "3 L^2 seed-H outputs plus 13 L^3 participating transfer outputs",
        "shared_final_decoder"=>channel == "x_only" ? "family.prefixes[end]; full-plaquette BP+OSD with declared p prior; syndrome-validated recovery" : "family.prefixes[end]; unit-weight data-edge PyMatching MWPM; syndrome-validated recovery",
        "failure"=>"either residual logical parity after final recovery",
        "joint_histogram_encoding"=>"mask = final_only + 2*full_prefix; each term is its state-failure Boolean",
        "point_rng"=>"Random.Xoshiro(point_seed)")
    for (key, value) in expected
        check!(counts, "protocol", get(protocol, key, nothing) == value, "protocol differs for $key")
    end
end

function check_sources!(counts, saved)
    required = ["Project.toml", "results/toric3d-h6-full-prefix-2026-10-15/run_study.jl"]
    isfile(joinpath(ROOT, "Manifest.toml")) && push!(required, "Manifest.toml")
    for (directory, _, files) in walkdir(joinpath(ROOT, "src")), file in files
        endswith(file, ".jl") && push!(required, relpath(joinpath(directory, file), ROOT))
    end
    check!(counts, "sources", Set(keys(saved)) == Set(required), "source manifest is incomplete or unexpected")
    observed = Dict{String,String}()
    for path in required
        observed[path] = digest(read(joinpath(ROOT, path)))
        check!(counts, "sources", observed[path] == saved[path], "source hash differs: $path")
    end
    return observed
end

function check_mirrors!(counts, directory, saved, bytes, snapshot)
    metadata_path, raw_path = joinpath(directory, "metadata.toml"), joinpath(directory, "raw.csv")
    metadata_bytes, raw_bytes = read(metadata_path), read(raw_path)
    metadata = TOML.parse(String(copy(metadata_bytes)))
    for key in ("format_version", "config", "runtime", "source_hashes", "protocol", "point_seeds_hex", "generation")
        check!(counts, "mirrors", metadata[key] == saved[key], "metadata differs for $key")
    end
    check!(counts, "mirrors", metadata["checkpoint_sha256"] == digest(bytes), "metadata checkpoint hash differs")
    check!(counts, "mirrors", metadata["raw_sha256"] == digest(raw_bytes), "metadata raw hash differs")
    check!(counts, "mirrors", metadata["current_snapshot"] == relpath(snapshot, directory), "metadata snapshot differs")
    lines = split(chomp(String(copy(raw_bytes))), '\n')
    fields = split(first(lines), ',')
    rows = [Dict(zip(fields, split(line, ','))) for line in lines[2:end]]
    grid = [(L, p) for L in saved["config"]["sizes"] for p in saved["config"]["rates"]]
    check!(counts, "mirrors", length(rows) == length(grid), "raw row count differs")
    for (index, (L, p)) in enumerate(grid)
        row, state = rows[index], saved["states"][index]
        check!(counts, "mirrors", row["channel"] == saved["config"]["channel"] &&
            all(parse(Int, row[name]) == L for name in ("steps", "Ly", "Lz")) &&
            parse(Float64, row["p"]) == p, "raw point identity differs")
        check!(counts, "mirrors", row["point_seed"] == saved["point_seeds_hex"][index], "raw point seed differs")
        check!(counts, "mirrors", parse(Int, row["shots"]) == state["completed_shots"], "raw completed shots differ")
        for (arm_index, arm) in enumerate(ARMS), (suffix, key) in
                (("logical_1_failures", "logical_1_counts"), ("logical_2_failures", "logical_2_counts"), ("failures", "failures"))
            check!(counts, "mirrors", parse(Int, row[arm*"_"*suffix]) == state[key][arm_index], "raw marginal differs")
        end
        for mask in 0:3
            check!(counts, "mirrors", parse(Int, row["joint_"*lpad(mask, 2, '0')]) == state["joint_counts"][mask+1],
                "raw joint bin differs")
        end
    end
    return Dict("metadata.toml"=>digest(metadata_bytes), "raw.csv"=>digest(raw_bytes))
end

function main(arguments=ARGS)
    length(arguments) == 4 || throw(ArgumentError(USAGE))
    started = time()
    directory = abspath(arguments[1])
    selected_L, selected_p = parse(Int, arguments[2]), parse(Float64, arguments[3])
    output = abspath(arguments[4])
    checkpoint_path = joinpath(directory, "checkpoint.toml")
    (output == directory || startswith(output, directory * "/")) &&
        throw(ArgumentError("verification output must be outside the read-only input directory"))
    counts = Dict{String,Int}()
    bytes = read(checkpoint_path)
    saved = TOML.parse(String(copy(bytes)))
    check!(counts, "checkpoint", saved["format_version"] == 1, "unsupported checkpoint format")
    payload = copy(saved)
    checksum = pop!(payload, "payload_sha256")
    check!(counts, "checkpoint", checksum == digest(canonical_toml(payload)), "payload checksum differs")
    snapshot_name = "checkpoint-" * lpad(saved["generation"], 6, '0') * "-" * checksum[1:12] * ".toml"
    snapshot = joinpath(directory, "checkpoints", snapshot_name)
    check!(counts, "checkpoint", read(snapshot) == bytes, "immutable checkpoint snapshot differs")
    config = saved["config"]
    channel = config["channel"]
    check!(counts, "checkpoint", channel in ("x_only", "z_only"), "unsupported noise channel")
    check!(counts, "checkpoint", config["dimensions"] == [[L, L, L] for L in config["sizes"]],
        "dimensions do not satisfy T=Ly=Lz=L")
    check_protocol!(counts, saved["protocol"], channel)
    mirror_hashes = check_mirrors!(counts, directory, saved, bytes, snapshot)
    grid = [(L, p) for L in config["sizes"] for p in config["rates"]]
    check!(counts, "checkpoint", length(grid) == length(saved["states"]), "grid and states differ")
    selected = findall(==((selected_L, selected_p)), grid)
    check!(counts, "checkpoint", length(selected) == 1, "requested point is absent or duplicated")
    index = only(selected)
    state = saved["states"][index]
    shots = state["completed_shots"]
    check!(counts, "checkpoint", 0 < shots <= config["shots"], "requested point has no saved samples")
    seeds = regenerate_seeds(config["seed"], length(grid))
    hex_seeds = [string(seed; base=16, pad=16) for seed in seeds]
    check!(counts, "seeds", all(occursin(r"^[0-9a-f]{16}$", seed) for seed in saved["point_seeds_hex"]),
        "point seeds are not 16 lowercase hexadecimal digits")
    check!(counts, "seeds", hex_seeds == saved["point_seeds_hex"], "regenerated point seed list differs")
    check!(counts, "seeds", parse(UInt64, saved["point_seeds_hex"][index]; base=16) == seeds[index], "chosen seed differs")
    observed_sources = check_sources!(counts, saved["source_hashes"])
    runtime = check_runtime!(counts, saved["runtime"])
    code = Toric3DEncodedSlice(selected_L, selected_L, selected_L)
    encoder = toric3d_encoded_encoder(code)
    check!(counts, "circuit", encoder.checkpoints == [7selected_L^2*r for r in 1:selected_L],
        "H+6 checkpoint schedule differs")
    check!(counts, "circuit", sum(length(operation.qubits) for operation in encoder.operations) == 13selected_L^3,
        "participating transfer-output count differs")
    check!(counts, "circuit", all(operation.gate in (:H, :CNOT) for operation in encoder.operations),
        "unexpected circuit operation")
    build_seconds = @elapsed begin
        if channel == "x_only"
            settings = Toric3DDecoderSettings()
            declared = Dict(string(k)=>v for (k,v) in pairs(TopoNoise._toric3d_settings_metadata(settings)))
            check!(counts, "decoder", declared == saved["protocol"]["decoder_settings"], "X decoder settings differ")
            decoder = build_toric3d_x_feedback_decoders(code; p=selected_p, settings)
        else
            check!(counts, "decoder", saved["protocol"]["decoder_settings"] ==
                Dict("weights"=>1.0, "merge_strategy"=>"smallest-weight"), "Z decoder settings differ")
            decoder = build_toric3d_z_feedback_decoders(code)
        end
    end
    geometry = independent_geometry(code, channel)
    joint, logical_1, logical_2, failures = zeros(Int, 4), zeros(Int, 2), zeros(Int, 2), zeros(Int, 2)
    rng = Xoshiro(seeds[index])
    replay_seconds = @elapsed for shot in 1:shots
        faults = regenerate_faults(rng, encoder, selected_p, channel)
        original = deepcopy(faults)
        mask = 0
        for (arm, feedback) in enumerate((false, true))
            # Direct public core calls. No runner import or second final decoder call.
            result = channel == "x_only" ?
                replay_toric3d_h6_full_prefix_x(encoder, faults, decoder; feedback, record_history=feedback) :
                replay_toric3d_h6_full_prefix_z(encoder, faults, decoder; feedback, record_history=feedback)
            check!(counts, "immutability", faults == original, "$(ARMS[arm]) mutated record $shot")
            check_final!(counts, result, geometry, channel)
            check_history!(counts, result, encoder, faults, code, geometry, channel; feedback)
            logical_1[arm] += result.logical_1
            logical_2[arm] += result.logical_2
            failures[arm] += result.failure
            result.failure && (mask |= 1 << (arm - 1))
        end
        joint[mask + 1] += 1
    end
    for mask in 0:3
        check!(counts, "counts", joint[mask + 1] == state["joint_counts"][mask + 1], "joint bin $mask differs")
    end
    for arm in 1:2
        check!(counts, "counts", logical_1[arm] == state["logical_1_counts"][arm], "logical-1 marginal differs")
        check!(counts, "counts", logical_2[arm] == state["logical_2_counts"][arm], "logical-2 marginal differs")
        check!(counts, "counts", failures[arm] == state["failures"][arm], "failure marginal differs")
        marginal = sum(joint[mask + 1] for mask in 0:3 if !iszero(mask & (1 << (arm - 1))))
        check!(counts, "counts", marginal == failures[arm], "joint marginal differs from failure count")
    end
    check!(counts, "counts", sum(joint) == shots, "joint histogram total differs")
    check!(counts, "read_only", read(checkpoint_path) == bytes, "input checkpoint changed during replay")
    check!(counts, "read_only", read(snapshot) == bytes, "input snapshot changed during replay")
    for (path, sha256) in mirror_hashes
        check!(counts, "read_only", digest(read(joinpath(directory, path))) == sha256, "input mirror changed during replay")
    end
    summary = Dict{String,Any}("status"=>"passed", "verification"=>"independent deterministic full-prefix point replay",
        "channel"=>channel, "dimensions"=>[selected_L, selected_L, selected_L], "p"=>selected_p,
        "shots"=>shots, "target_shots"=>config["shots"], "point_index"=>index, "base_seed"=>config["seed"],
        "point_seed_hex"=>hex_seeds[index], "regenerated_point_seeds_hex"=>hex_seeds, "arm_order"=>collect(ARMS),
        "joint_counts"=>joint, "logical_1_counts"=>logical_1, "logical_2_counts"=>logical_2, "failures"=>failures,
        "assertions_total"=>sum(values(counts)), "assertions_by_category"=>counts,
        "histories_checked"=>shots, "prefixes_checked"=>shots*selected_L, "last_prefixes_checked"=>shots,
        "final_outcomes_checked"=>2shots, "build_seconds"=>build_seconds, "replay_seconds"=>replay_seconds,
        "total_seconds"=>time()-started, "input_directory"=>directory, "input_checkpoint"=>checkpoint_path,
        "input_checkpoint_sha256"=>digest(bytes), "input_payload_sha256"=>checksum, "input_snapshot"=>snapshot,
        "input_generation"=>saved["generation"], "input_mirror_sha256"=>mirror_hashes, "source_hashes"=>observed_sources,
        "verifier_sha256"=>digest(read(@__FILE__)), "runtime"=>runtime,
        "final_validation"=>"independent physical-check and logical-support matrix products; no extra final decode",
        "fault_validation"=>"independent scalar Xoshiro draws; X includes its established p_z=0 draws; Z draws seed first",
        "history_validation"=>"independent H/CNOT propagation and recorded physical recoveries; all physical prefix checks; Q0:Qr support and no future correction at every r=1:T; no separate Z seed correction",
        "independence_limit"=>"shares public circuit/check/logical-support definitions and decoder implementations; does not import runner, production sampler, propagation helper or prefix-map builder")
    mkpath(dirname(output))
    open(output * ".tmp", "w") do io
        TOML.print(io, summary; sorted=true)
    end
    mv(output * ".tmp", output; force=true)
    println("Passed ", channel, " L=", selected_L, " p=", selected_p, " shots=", shots,
        ": ", summary["assertions_total"], " assertions; joint=", joint)
    println("Saved ", output)
    return summary
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    H6FullPrefixIndependentReplay.main()
end
