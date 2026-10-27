"""The literal sequential Clifford schedule and full-prefix checkpoints."""
struct Toric3DEncoder
    code::ToricCode3D
    operations::Vector{EncoderOperation}
    checkpoints::Vector{Int}
    prefixes::Vector{ToricCode3D}
    prefix_qubits::Vector{Vector{Int}}
end

"""Prepare `|0_L⟩` by visiting vertices in x, y, z order."""
function toric3d_encoder(code::ToricCode3D)
    operations = EncoderOperation[]
    checkpoints = Int[]
    prefixes = ToricCode3D[]
    prefix_qubits = Vector{Int}[]
    q(d, x, y, z) = data_qubit_index(code, d, x, mod(y, code.Ly), mod(z, code.Lz))
    for x in 1:code.Lx
        for y in 0:code.Ly-1, z in 0:code.Lz-1
            control = q(:x, x, y, z)
            push!(operations, EncoderOperation(:H, [control]))
            for target in (q(:x, x - 1, y, z),
                           q(:y, x, y, z), q(:y, x, y - 1, z),
                           q(:z, x, y, z), q(:z, x, y, z - 1))
                push!(operations, EncoderOperation(:CNOT, [control, target]))
            end
        end
        push!(checkpoints, length(operations))
        prefix = toric3d_prefix(code, x)
        push!(prefixes, prefix.code)
        push!(prefix_qubits, prefix.qubits)
    end
    return Toric3DEncoder(code, operations, checkpoints, prefixes, prefix_qubits)
end

function _toric3d_pulled_to_zero_input(
        n::Int, operations::AbstractVector{EncoderOperation},
        pauli::Symbol, support::AbstractVector{Int})
    x, z = falses(n), falses(n)
    negative = Ref(false)
    if pauli === :X
        x[support] .= true
    else
        z[support] .= true
    end
    for operation in Iterators.reverse(operations)
        operation.gate in (:H, :CNOT) || return false
        _conjugate_pauli!(x, z, negative, operation)
    end
    return !any(x) && !negative[]
end

"""Algebraically check every prefix and final prepared stabilizer/logical Z."""
function verify_toric3d_encoder(encoder::Toric3DEncoder)
    code = encoder.code
    length(encoder.checkpoints) == code.Lx || return false
    length(encoder.prefixes) == length(encoder.prefix_qubits) == code.Lx ||
        return false
    n = data_qubit_count(code)
    for j in 1:code.Lx
        checkpoint = encoder.checkpoints[j]
        1 <= checkpoint <= length(encoder.operations) || return false
        expected = toric3d_prefix(code, j)
        prefix = encoder.prefixes[j]
        mapping = encoder.prefix_qubits[j]
        (prefix.Lx, prefix.Ly, prefix.Lz) ==
            (expected.code.Lx, expected.code.Ly, expected.code.Lz) || return false
        mapping == expected.qubits || return false
        operations = @view encoder.operations[1:checkpoint]
        all(all(q in mapping for q in operation.qubits) for operation in operations) ||
            return false
        for (pauli, checks) in ((:X, star_checks(prefix)),
                                (:Z, plaquette_checks(prefix)),
                                (:Z, [logical_z_support(prefix)]))
            all(_toric3d_pulled_to_zero_input(n, operations, pauli, mapping[support])
                for support in checks) || return false
        end
    end
    return encoder.checkpoints[end] == length(encoder.operations)
end
