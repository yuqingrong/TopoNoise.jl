"""Independent post-gate X/Z bits on each active output of a fixed schedule."""
struct Toric3DFaultRecord
    Lx::Int
    Ly::Int
    Lz::Int
    x::BitMatrix
    z::BitMatrix

    function Toric3DFaultRecord(Lx::Integer, Ly::Integer, Lz::Integer,
                                x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool})
        ToricCode3D(Lx, Ly, Lz)
        size(x) == size(z) && size(x, 1) == 2 ||
            throw(ArgumentError("fault matrices must have matching 2×operation shapes"))
        return new(Int(Lx), Int(Ly), Int(Lz), BitMatrix(x), BitMatrix(z))
    end
end

Base.:(==)(a::Toric3DFaultRecord, b::Toric3DFaultRecord) =
    (a.Lx, a.Ly, a.Lz) == (b.Lx, b.Ly, b.Lz) && a.x == b.x && a.z == b.z

function _validate_toric3d_fault_record(encoder::Toric3DEncoder,
                                        record::Toric3DFaultRecord)
    code = encoder.code
    (record.Lx, record.Ly, record.Lz) == (code.Lx, code.Ly, code.Lz) ||
        throw(ArgumentError("fault record has incompatible code dimensions"))
    size(record.x) == size(record.z) == (2, length(encoder.operations)) ||
        throw(ArgumentError("fault record has incompatible operation count"))
    for (step, operation) in enumerate(encoder.operations)
        if operation.gate === :H
            !(record.x[2, step] || record.z[2, step]) ||
                throw(ArgumentError("H gates have no second active output"))
        elseif operation.gate !== :CNOT
            throw(ArgumentError("unsupported encoder gate"))
        end
    end
    return nothing
end

"""Sample X and Z Bernoulli faults after every H/CNOT active output."""
function sample_toric3d_faults(rng::Random.AbstractRNG,
                               encoder::Toric3DEncoder;
                               p_x::Real=0.0, p_z::Real=0.0)
    _validate_fault_rate(p_x, "p_x")
    _validate_fault_rate(p_z, "p_z")
    count = length(encoder.operations)
    x, z = falses(2, count), falses(2, count)
    for (step, operation) in enumerate(encoder.operations)
        for slot in 1:length(operation.qubits)
            x[slot, step] = rand(rng) < p_x
            z[slot, step] = rand(rng) < p_z
        end
    end
    return Toric3DFaultRecord(encoder.code.Lx, encoder.code.Ly,
                              encoder.code.Lz, x, z)
end

"""Return a new record with a composed `:X`, `:Z`, or `:Y` output fault."""
function with_toric3d_fault(encoder::Toric3DEncoder,
                            record::Toric3DFaultRecord,
                            step::Integer, qubit::Integer, pauli::Symbol)
    _validate_toric3d_fault_record(encoder, record)
    pauli in (:X, :Z, :Y) ||
        throw(ArgumentError("pauli must be :X, :Z, or :Y"))
    1 <= step <= length(encoder.operations) ||
        throw(ArgumentError("fault step is out of range"))
    slot = findfirst(==(Int(qubit)), encoder.operations[Int(step)].qubits)
    isnothing(slot) && throw(ArgumentError("qubit is not active at this gate"))
    x, z = copy(record.x), copy(record.z)
    pauli in (:X, :Y) && (x[slot, step] = !x[slot, step])
    pauli in (:Z, :Y) && (z[slot, step] = !z[slot, step])
    return Toric3DFaultRecord(record.Lx, record.Ly, record.Lz, x, z)
end

"""Advance a mutable frame through operations `first:last` and their faults."""
function propagate_toric3d_segment!(frame::PauliFrame,
                                    encoder::Toric3DEncoder,
                                    record::Toric3DFaultRecord,
                                    first::Integer, last::Integer)
    _validate_toric3d_fault_record(encoder, record)
    n = data_qubit_count(encoder.code)
    length(frame.x) == n || throw(ArgumentError("Pauli frame has the wrong size"))
    count = length(encoder.operations)
    1 <= first <= count + 1 && 0 <= last <= count && first <= last + 1 ||
        throw(ArgumentError("operation segment is out of range"))
    for step in first:last
        operation = encoder.operations[step]
        _propagate_frame_operation!(frame.x, frame.z, operation)
        for (slot, qubit) in enumerate(operation.qubits)
            frame.x[qubit] = xor(frame.x[qubit], record.x[slot, step])
            frame.z[qubit] = xor(frame.z[qubit], record.z[slot, step])
        end
    end
    return frame
end

"""Measure X stars against Z faults and Z plaquettes against X faults."""
function toric3d_syndrome(code::ToricCode3D, frame::PauliFrame)
    length(frame.x) == data_qubit_count(code) ||
        throw(ArgumentError("Pauli frame has the wrong size"))
    return (
        stars=BitVector(isodd(count(frame.z[s])) for s in star_checks(code)),
        plaquettes=BitVector(isodd(count(frame.x[s])) for s in plaquette_checks(code)),
    )
end

"""Both logical Pauli errors, named by the error operator."""
function toric3d_logical_parities(code::ToricCode3D, frame::PauliFrame)
    length(frame.x) == data_qubit_count(code) ||
        throw(ArgumentError("Pauli frame has the wrong size"))
    return (
        logical_x=isodd(count(frame.x[logical_z_support(code)])),
        logical_z=isodd(count(frame.z[logical_x_support(code)])),
    )
end
