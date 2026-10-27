"""Fixed Clifford transfer schedule with one checkpoint per completed slice."""
struct Toric3DSliceEncoder
    code::Toric3DSliceLattice
    operations::Vector{EncoderOperation}
    checkpoints::Vector{Int}
end

"""
Build the diagonal transfer circuit. A fresh GHZ triad supplies one vertex
amplitude; old-slice pair parities K12 and K13 control only new targets. Prepared
qubits are never Hadamard targets or CNOT targets in a later slice.
"""
function toric3d_slice_encoder(code::Toric3DSliceLattice)
    operations = EncoderOperation[]
    checkpoints = Int[]
    for t in 1:code.steps
        for u in 0:code.Ly-1, v in 0:code.Lz-1
            new = [data_qubit_index(code, i, t, u, v) for i in 1:3]
            push!(operations, EncoderOperation(:H, [new[1]]),
                EncoderOperation(:CNOT, [new[1], new[2]]),
                EncoderOperation(:CNOT, [new[1], new[3]]))
            for j in (2, 3)
                for control in _toric3d_slice_pair(code, t, u, v, 1, j)
                    push!(operations, EncoderOperation(:CNOT, [control, new[j]]))
                end
            end
        end
        push!(checkpoints, length(operations))
    end
    return Toric3DSliceEncoder(code, operations, checkpoints)
end

"""
Verify all prefix stars, physical plaquettes, initial boundary checks and
current-slice compatibility checks by signed binary Clifford conjugation back
to the all-zero input. This includes the final three-edge X frontier stars.
"""
function verify_toric3d_slice_encoder(encoder::Toric3DSliceEncoder)
    code = encoder.code
    cells = code.Ly * code.Lz
    encoder.checkpoints == [7 * cells * t for t in 1:code.steps] || return false
    length(encoder.operations) == encoder.checkpoints[end] || return false
    n = data_qubit_count(code)
    for operation in encoder.operations
        required = operation.gate === :H ? 1 : operation.gate === :CNOT ? 2 : 0
        required > 0 && length(operation.qubits) == required || return false
        all(1 <= q <= n for q in operation.qubits) || return false
        length(unique(operation.qubits)) == required || return false
    end
    for t in 1:code.steps
        prefix = Toric3DSliceLattice(t, code.Ly, code.Lz)
        operations = @view encoder.operations[1:encoder.checkpoints[t]]
        all(all(q <= data_qubit_count(prefix) for q in operation.qubits)
            for operation in operations) || return false
        for (pauli, checks) in ((:X, star_checks(prefix)),
                                (:Z, plaquette_checks(prefix)),
                                (:Z, toric3d_slice_boundary_checks(prefix)),
                                (:Z, toric3d_slice_checks(prefix, t)))
            all(_toric3d_pulled_to_zero_input(n, operations, pauli, support)
                for support in checks) || return false
        end
    end
    return true
end
