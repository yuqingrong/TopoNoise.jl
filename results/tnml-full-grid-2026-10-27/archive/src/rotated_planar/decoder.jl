"""Stable Julia-owned handle to the two CSS PyMatching decoders."""
struct MatchingDecoders
    x_error_matcher::PythonCall.Py
    z_error_matcher::PythonCall.Py
    backend::Symbol
    backend_version::VersionNumber
    a_s_detector_count::Int
    b_p_detector_count::Int
end

function _validate_graphlike_check_matrix(matrix, family::Symbol)
    column_weights = vec(sum(matrix; dims=1))
    all(<=(2), column_weights) ||
        throw(ArgumentError(
            "$family check matrix is not graphlike: every column must have weight at most two"))
    return nothing
end

function _load_pymatching()
    try
        return PythonCall.pyimport("pymatching")
    catch error
        message = sprint(showerror, error)
        throw(ErrorException(
            "PyMatching could not be imported. Run `import CondaPkg; CondaPkg.resolve()` " *
            "for this project, then retry. Original error: $message"))
    end
end

function _numpy_binary_matrix(numpy, matrix::AbstractMatrix{Bool})
    return numpy.asarray(Matrix{UInt8}(matrix); dtype=numpy.uint8)
end

function _fault_observable(qubit_count::Int, support::AbstractVector{<:Integer})
    observable = falses(1, qubit_count)
    observable[1, support] .= true
    return observable
end

"""
    build_matching_decoders(code) -> MatchingDecoders

Build independent final-syndrome PyMatching decoders for X and Z errors.
The A_s/Z-check sector decodes X faults against logical Z, while the
B_p/X-check sector decodes Z faults against logical X.
"""
function build_matching_decoders(code::RotatedPlanarCode)
    a_s = a_s_check_matrix(code)
    b_p = b_p_check_matrix(code)
    _validate_graphlike_check_matrix(a_s, :A_s)
    _validate_graphlike_check_matrix(b_p, :B_p)

    pymatching = _load_pymatching()
    numpy = PythonCall.pyimport("numpy")
    qubit_count = data_qubit_count(code)
    x_error_matcher = pymatching.Matching.from_check_matrix(
        _numpy_binary_matrix(numpy, a_s);
        weights=1.0,
        faults_matrix=_numpy_binary_matrix(
            numpy, _fault_observable(qubit_count, logical_z_support(code))),
    )
    z_error_matcher = pymatching.Matching.from_check_matrix(
        _numpy_binary_matrix(numpy, b_p);
        weights=1.0,
        faults_matrix=_numpy_binary_matrix(
            numpy, _fault_observable(qubit_count, logical_x_support(code))),
    )
    backend_version = VersionNumber(PythonCall.pyconvert(String, pymatching.__version__))
    return MatchingDecoders(
        x_error_matcher, z_error_matcher, :PyMatching, backend_version,
        size(a_s, 1), size(b_p, 1))
end

function _prediction_bits(predictions)
    numpy = PythonCall.pyimport("numpy")
    values = numpy.ascontiguousarray(predictions; dtype=numpy.uint8).reshape(-1)
    bytes = PythonCall.pyconvert(Vector{UInt8}, values.tobytes())
    return BitVector(bytes .!= 0)
end

"""
    decode_logical_parities(decoders, a_s_syndromes, b_p_syndromes)

Decode a batch using only the final A_s and B_p syndrome matrices. Rows are
shots and columns are checks. The returned vectors are predicted logical-X
and logical-Z parities, respectively; actual faults never cross this boundary.
"""
function decode_logical_parities(
        decoders::MatchingDecoders,
        a_s_syndromes::AbstractMatrix{Bool},
        b_p_syndromes::AbstractMatrix{Bool})
    size(a_s_syndromes, 2) == decoders.a_s_detector_count ||
        throw(DimensionMismatch(
            "A_s syndrome matrix must have $(decoders.a_s_detector_count) columns"))
    size(b_p_syndromes, 2) == decoders.b_p_detector_count ||
        throw(DimensionMismatch(
            "B_p syndrome matrix must have $(decoders.b_p_detector_count) columns"))
    size(a_s_syndromes, 1) == size(b_p_syndromes, 1) ||
        throw(DimensionMismatch("A_s and B_p syndrome matrices must have equal row counts"))

    numpy = PythonCall.pyimport("numpy")
    predicted_x = decoders.x_error_matcher.decode_batch(
        _numpy_binary_matrix(numpy, a_s_syndromes);
        enable_correlations=false,
    )
    predicted_z = decoders.z_error_matcher.decode_batch(
        _numpy_binary_matrix(numpy, b_p_syndromes);
        enable_correlations=false,
    )
    return (
        logical_x=_prediction_bits(predicted_x),
        logical_z=_prediction_bits(predicted_z),
    )
end

"""One circuit-level Monte Carlo estimate with code, noise, and uncertainty metadata."""
struct LogicalFailurePoint
    distance::Int
    boundary_orientation::Symbol
    construction::Symbol
    logical_state::Symbol
    clock::Symbol
    p_x::Float64
    p_z::Float64
    shots::Int
    seed::Union{Nothing,Int}
    logical_x_failures::Int
    logical_x_failure_rate::Float64
    logical_x_standard_error::Float64
    logical_z_failures::Int
    logical_z_failure_rate::Float64
    logical_z_standard_error::Float64
    any_logical_failures::Int
    any_logical_failure_rate::Float64
    any_logical_standard_error::Float64
    state_failures::Int
    state_failure_rate::Float64
    state_standard_error::Float64
end

function _positive_machine_int(value::Integer, name::AbstractString)
    value > 0 || throw(ArgumentError("$name must be positive"))
    try
        return Int(value)
    catch
        throw(ArgumentError("$name must be representable as Int"))
    end
end

function _seed_metadata(seed)
    seed === nothing && return nothing
    try
        return Int(seed)
    catch
        throw(ArgumentError("seed must be nothing or convertible to Int"))
    end
end

function _propagate_frame_batch_operation!(
        x::BitMatrix, z::BitMatrix, operation::EncoderOperation)
    if operation.gate === :H
        qubit = only(operation.qubits)
        previous_x = copy(@view x[:, qubit])
        @views x[:, qubit] .= z[:, qubit]
        @views z[:, qubit] .= previous_x
    elseif operation.gate === :X
        nothing
    elseif operation.gate === :Z
        nothing
    elseif operation.gate === :CNOT
        control_qubit, target_qubit = operation.qubits
        @views x[:, target_qubit] .= xor.(
            x[:, target_qubit], x[:, control_qubit])
        @views z[:, control_qubit] .= xor.(
            z[:, control_qubit], z[:, target_qubit])
    else
        error("unsupported encoder gate $(operation.gate)")
    end
    return nothing
end

function _inject_random_frame_batch!(
        rng::Random.AbstractRNG, x::BitMatrix, z::BitMatrix,
        support::AbstractVector{<:Integer}, noise::CircuitPauliNoise)
    shot_count = size(x, 1)
    x_events = rand(rng, shot_count, length(support)) .< noise.p_x
    z_events = rand(rng, shot_count, length(support)) .< noise.p_z
    @views x[:, support] .= xor.(x[:, support], x_events)
    @views z[:, support] .= xor.(z[:, support], z_events)
    return nothing
end

function _sample_frame_batch(
        rng::Random.AbstractRNG, encoder::PlaquetteEncoder,
        noise::CircuitPauliNoise, shot_count::Int)
    qubit_count = data_qubit_count(encoder.code)
    x = falses(shot_count, qubit_count)
    z = falses(shot_count, qubit_count)

    if noise.clock === :post_encoding
        # The reference state is the ideally encoded state. Draw one fault
        # per data qubit, with no subsequent Clifford propagation.
        _inject_random_frame_batch!(rng, x, z, collect(1:qubit_count), noise)
    elseif noise.clock === :gate_layer
        for layer in gate_layers(encoder)
            for operation in layer.operations
                _propagate_frame_batch_operation!(x, z, operation)
            end
            _inject_random_frame_batch!(
                rng, x, z, layer.active_qubits, noise)
        end
    else
        for block in plaquette_blocks(encoder)
            for layer in block.layers, operation in layer.operations
                _propagate_frame_batch_operation!(x, z, operation)
            end
            block.kind === :plaquette && _inject_random_frame_batch!(
                rng, x, z, block.source_support, noise)
        end
    end
    return x, z
end

function _support_parities(frames::BitMatrix, supports)
    parities = falses(size(frames, 1), length(supports))
    for (column, support) in enumerate(supports)
        @views parities[:, column] .= isodd.(vec(sum(frames[:, support]; dims=2)))
    end
    return parities
end

function _logical_parities(frames::BitMatrix, support::AbstractVector{<:Integer})
    @views return BitVector(isodd.(vec(sum(frames[:, support]; dims=2))))
end

_binomial_standard_error(rate::Float64, shots::Int) =
    sqrt(rate * (1 - rate) / shots)

"""
    estimate_logical_failure(rng, encoder, noise; shots=10_000,
                             batch_size=10_000, seed=nothing)

Estimate logical failure probabilities with chunked shot-by-qubit binary
Pauli frames. The supplied RNG is the only source of randomness; `seed` is
retained solely as result metadata.
"""
function estimate_logical_failure(
        rng::Random.AbstractRNG,
        encoder::PlaquetteEncoder,
        noise::CircuitPauliNoise;
        shots::Integer=10_000,
        batch_size::Integer=10_000,
        seed=nothing,
    )::LogicalFailurePoint
    shot_count = _positive_machine_int(shots, "shots")
    chunk_limit = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)

    code = encoder.code
    decoders = build_matching_decoders(code)
    logical_x_failures = 0
    logical_z_failures = 0
    any_logical_failures = 0
    remaining = shot_count

    while remaining > 0
        chunk_size = min(remaining, chunk_limit)
        x_frames, z_frames = _sample_frame_batch(
            rng, encoder, noise, chunk_size)
        a_s_syndromes = _support_parities(x_frames, a_s_checks(code))
        b_p_syndromes = _support_parities(z_frames, b_p_checks(code))
        actual_x = _logical_parities(x_frames, logical_z_support(code))
        actual_z = _logical_parities(z_frames, logical_x_support(code))

        # Only final syndrome matrices cross into Python. Actual logical
        # parities remain in Julia and are compared after batch decoding.
        predicted = decode_logical_parities(
            decoders, a_s_syndromes, b_p_syndromes)
        failed_x = xor.(predicted.logical_x, actual_x)
        failed_z = xor.(predicted.logical_z, actual_z)
        logical_x_failures += count(failed_x)
        logical_z_failures += count(failed_z)
        any_logical_failures += count(failed_x .| failed_z)
        remaining -= chunk_size
    end

    state_failures = encoder.logical_state in (:zero, :one) ?
        logical_x_failures : logical_z_failures
    logical_x_rate = logical_x_failures / shot_count
    logical_z_rate = logical_z_failures / shot_count
    any_logical_rate = any_logical_failures / shot_count
    state_rate = state_failures / shot_count
    return LogicalFailurePoint(
        distance(code), boundary_orientation(code), encoder.construction,
        encoder.logical_state, noise.clock, noise.p_x, noise.p_z,
        shot_count, seed_value,
        logical_x_failures, logical_x_rate,
        _binomial_standard_error(logical_x_rate, shot_count),
        logical_z_failures, logical_z_rate,
        _binomial_standard_error(logical_z_rate, shot_count),
        any_logical_failures, any_logical_rate,
        _binomial_standard_error(any_logical_rate, shot_count),
        state_failures, state_rate,
        _binomial_standard_error(state_rate, shot_count),
    )
end
