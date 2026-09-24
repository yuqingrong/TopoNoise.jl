# Export the current Bp schedules and independent small-distance Yao oracles.
using TopoNoise, TOML, SHA, Random, LinearAlgebra
import Yao

const RUN = @__DIR__
const ROOT = normpath(joinpath(RUN, "..", ".."))

function write_once(path, value)
    text = sprint(io -> TOML.print(io, value; sorted=true))
    if isfile(path)
        read(path, String) == text || error("Refusing to replace a different archived export: $path")
    else
        write(path, text)
    end
end

function export_schedule(d)
    code = RotatedPlanarCode(d; boundary_orientation=:x_ns)
    enc = rotated_planar_encoder(code; construction=:bp)
    @assert verify_encoder_tableau(enc)
    @assert enc.input_state === :zero && enc.logical_state === :zero
    return Dict("distance" => d, "construction" => "bp", "boundary_orientation" => "x_ns",
        "input_state" => "zero", "clock" => "gate_layer", "a_s" => a_s_checks(code),
        "b_p" => b_p_checks(code), "logical_x" => logical_x_support(code),
        "logical_z" => logical_z_support(code),
        "layers" => [Dict("active_qubits" => layer.active_qubits,
            "operations" => [Dict("gate" => String(op.gate), "qubits" => op.qubits)
                             for op in layer.operations]) for layer in gate_layers(enc)])
end

function oracle_fixtures()
    code = RotatedPlanarCode(3)
    enc = rotated_planar_encoder(code; construction=:bp)
    layers = gate_layers(enc)
    locations = [(i, q) for (i, layer) in enumerate(layers) for q in layer.active_qubits]
    n = data_qubit_count(code) + 1
    initial = Yao.zero_state(n)
    for layer in layers, op in layer.operations
        Yao.apply!(initial, TopoNoise._yao_operation(n, op))
    end
    Yao.apply!(initial, Yao.put(n, n => Yao.H))
    for q in logical_x_support(code)
        Yao.apply!(initial, Yao.control(n, n, q => Yao.X))
    end
    for layer in reverse(layers), op in reverse(layer.operations)
        Yao.apply!(initial, TopoNoise._yao_operation(n, op))
    end
    decoders = build_matching_decoders(code)
    fixtures = Any[]
    rng = MersenneTwister(20260923)
    for (channel, gate) in (("y_only", "Y"), ("hadamard", "H"))
        patterns = [[i] for i in eachindex(locations)]
        append!(patterns, [sort(randperm(rng, length(locations))[1:rand(rng, 2:6)]) for _ in 1:16])
        gate == "H" && push!(patterns, [5, 5])
        for (i, indices) in enumerate(patterns)
            faults = [Dict("step" => locations[k][1], "qubit" => locations[k][2], "gate" => gate)
                      for k in indices]
            reg = copy(initial)
            for (step, layer) in enumerate(layers)
                for op in layer.operations
                    Yao.apply!(reg, TopoNoise._yao_operation(n, op))
                end
                for fault in faults
                    fault["step"] == step || continue
                    Yao.apply!(reg, Yao.put(n, fault["qubit"] => (gate == "Y" ? Yao.Y : Yao.H)))
                end
            end
            state = Yao.statevec(reg)
            @assert isapprox(norm(state), 1; atol=1e-12)
            fixture = Dict{String,Any}("name" => "$channel-$i", "channel" => channel,
                "faults" => faults, "state_real" => real.(state), "state_imag" => imag.(state))
            if gate == "Y"
                record = sample_fault_record(MersenneTwister(0), enc, CircuitPauliNoise(0))
                for fault in faults
                    record = with_pauli_fault(record, fault["step"], fault["qubit"], :Y)
                end
                frame = propagate_pauli_frame(enc, record)
                syndrome = measure_syndrome(code, frame)
                lx = isodd(count(frame.x[logical_z_support(code)]))
                lz = isodd(count(frame.z[logical_x_support(code)]))
                predicted = decode_logical_parities(decoders,
                    reshape(syndrome.a_s, 1, :), reshape(syndrome.b_p, 1, :))
                fixture["raw_measurements"] = vcat(syndrome.a_s, syndrome.b_p, [lx, lz])
                fixture["logical_failures"] = [xor(lx, only(predicted.logical_x)),
                                                xor(lz, only(predicted.logical_z))]
            end
            push!(fixtures, fixture)
        end
    end
    return Dict("fixtures" => fixtures, "oracle" => "Yao statevector and TopoNoise Pauli propagation")
end

function snapshot_sources()
    paths = ["Project.toml", "Manifest.toml", "CondaPkg.toml"]
    for (directory, _, files) in walkdir(joinpath(ROOT, "src")), file in files
        endswith(file, ".jl") && push!(paths, relpath(joinpath(directory, file), ROOT))
    end
    hashes = Dict{String,String}()
    for relative in sort(paths)
        source = joinpath(ROOT, relative)
        isfile(source) || continue
        data = read(source)
        hashes[relative] = bytes2hex(sha256(data))
        target = joinpath(RUN, "source_snapshot", relative)
        mkpath(dirname(target))
        if isfile(target)
            read(target) == data || error("Source changed since archive creation: $relative")
        else
            write(target, data)
        end
    end
    write_once(joinpath(RUN, "source_hashes.toml"), Dict("sha256" => hashes,
        "julia_version" => string(VERSION)))
end

snapshot_sources()
write_once(joinpath(RUN, "circuits.toml"), Dict("circuits" => export_schedule.([3, 9, 11, 13, 15])))
write_once(joinpath(RUN, "oracle_fixtures.toml"), oracle_fixtures())
println("Exported five verified Bp schedules and 73 independent Yao fixtures.")
