"""A validated unit-weight MWPM decoder for one local Z-feedback graph."""
struct Toric3DZFeedbackDecoder
    window::Toric3DZFeedbackWindow
    matrix::SparseArrays.SparseMatrixCSC{UInt8,Int}
    matcher::PythonCall.Py
    version::String
end

"""Cached seed, transfer-window, and growing-prefix Z decoders."""
struct Toric3DZFeedbackDecoderFamily
    code::Toric3DEncodedSlice
    seed::Toric3DZFeedbackDecoder
    transfer::Toric3DZFeedbackDecoder
    prefixes::Vector{Toric3DZThresholdDecoder}
end

function _build_toric3d_z_feedback_decoder(window::Toric3DZFeedbackWindow,
        pymatching, version::String)
    matrix = SparseArrays.sparse(UInt8.(window.matrix))
    all(==(2),vec(sum(matrix;dims=1))) ||
        error("Z-feedback recovery graph must give every edge two endpoints")
    n = size(matrix,2)
    matcher = pymatching.Matching.from_check_matrix(_toric3d_scipy_sparse(matrix);
        weights=1.0,merge_strategy="smallest-weight",min_num_fault_ids=n)
    return Toric3DZFeedbackDecoder(window,matrix,matcher,version)
end

function build_toric3d_z_feedback_decoders(code::Toric3DEncodedSlice)
    pymatching = _load_pymatching()
    version = _toric3d_module_version(pymatching,"pymatching")
    seed = _build_toric3d_z_feedback_decoder(
        toric3d_z_feedback_window(code,0),pymatching,version)
    transfer_window = toric3d_z_feedback_window(code,1)
    for r in 2:code.lattice.steps
        toric3d_z_feedback_window(code,r).matrix == transfer_window.matrix ||
            error("ordered Z-feedback transfer windows are not reusable")
    end
    transfer = _build_toric3d_z_feedback_decoder(transfer_window,pymatching,version)
    prefixes = [build_toric3d_z_threshold_decoder(
        toric3d_encoded_prefix(code,t).code) for t in 1:code.lattice.steps]
    return Toric3DZFeedbackDecoderFamily(code,seed,transfer,prefixes)
end

function _toric3d_validate_z_feedback_recovery(decoder::Toric3DZFeedbackDecoder,
        syndrome::AbstractVector{Bool}, recovery::AbstractVector{Bool})
    length(syndrome) == size(decoder.matrix,1) ||
        throw(DimensionMismatch("Z-feedback syndrome has the wrong length"))
    length(recovery) == size(decoder.matrix,2) ||
        throw(DimensionMismatch("Z-feedback recovery has the wrong length"))
    checked = BitVector(recovery)
    BitVector(isodd.(decoder.matrix*Vector{UInt8}(checked))) == BitVector(syndrome) ||
        error("invalid Z-feedback recovery: its syndrome differs from the measurement")
    return checked
end

function _toric3d_validate_z_feedback_batch(decoder::Toric3DZFeedbackDecoder,
        syndromes::AbstractMatrix{Bool}, recoveries::AbstractMatrix{Bool})
    size(recoveries) == (size(decoder.matrix,2),size(syndromes,2)) ||
        throw(DimensionMismatch("Z-feedback recovery batch has the wrong dimensions"))
    BitMatrix(isodd.(decoder.matrix*Matrix{UInt8}(recoveries))) == syndromes ||
        error("invalid Z-feedback recovery batch: its syndrome differs from the measurement")
    return BitMatrix(recoveries)
end

"""Decode one local measured syndrome and reject invalid graph recovery."""
function decode_toric3d_z_feedback(decoder::Toric3DZFeedbackDecoder,
        syndrome::AbstractVector{Bool})::BitVector
    length(syndrome) == size(decoder.matrix,1) ||
        throw(DimensionMismatch("Z-feedback syndrome has the wrong length"))
    iseven(count(syndrome)) ||
        throw(ArgumentError("closed Z-feedback graph syndromes must have even parity"))
    n = size(decoder.matrix,2)
    !any(syndrome) && return falses(n)
    numpy = PythonCall.pyimport("numpy")
    predicted = decoder.matcher.decode(_toric3d_numpy_vector(numpy,syndrome))
    recovery = _toric3d_py_bits(numpy,predicted,n)
    return _toric3d_validate_z_feedback_recovery(decoder,syndrome,recovery)
end

"""Decode a columns-as-shots syndrome batch and validate every recovery."""
function decode_toric3d_z_feedback_batch(decoder::Toric3DZFeedbackDecoder,
        syndromes::AbstractMatrix{Bool})::BitMatrix
    size(syndromes,1) == size(decoder.matrix,1) ||
        throw(DimensionMismatch("Z-feedback syndrome batch has the wrong dimensions"))
    all(iseven,vec(sum(syndromes;dims=1))) ||
        throw(ArgumentError("closed Z-feedback graph syndromes must have even parity"))
    n,shots = size(decoder.matrix,2),size(syndromes,2)
    shots == 0 && return falses(n,0)
    numpy = PythonCall.pyimport("numpy")
    input = numpy.asarray(Matrix{UInt8}(permutedims(syndromes));dtype=numpy.uint8)
    predicted = decoder.matcher.decode_batch(input)
    recovery = BitMatrix(reshape(_toric3d_py_bits(numpy,predicted,n*shots),n,shots))
    return _toric3d_validate_z_feedback_batch(decoder,syndromes,recovery)
end
