"""Reproducible BP+OSD controls for X recovery in the 3D slab."""
struct Toric3DDecoderSettings
    max_iter::Int
    bp_method::String
    ms_scaling_factor::Float64
    schedule::String
    omp_thread_count::Int
    osd_method::String
    osd_order::Int

    function Toric3DDecoderSettings(; max_iter::Integer=50,
            bp_method::AbstractString="minimum_sum",
            ms_scaling_factor::Real=1.0,
            schedule::AbstractString="parallel",
            omp_thread_count::Integer=1,
            osd_method::AbstractString="OSD_CS",
            osd_order::Integer=2)
        0 <= max_iter <= typemax(Int) ||
            throw(ArgumentError("max_iter must be a nonnegative Int"))
        bp_method in ("minimum_sum", "product_sum") ||
            throw(ArgumentError("bp_method must be minimum_sum or product_sum"))
        isfinite(ms_scaling_factor) && ms_scaling_factor >= 0 ||
            throw(ArgumentError("ms_scaling_factor must be finite and nonnegative"))
        schedule in ("parallel", "serial") ||
            throw(ArgumentError("schedule must be parallel or serial"))
        1 <= omp_thread_count <= typemax(Int) ||
            throw(ArgumentError("omp_thread_count must be a positive Int"))
        osd_method in ("OSD_0", "OSD_E", "OSD_CS") ||
            throw(ArgumentError("osd_method must be OSD_0, OSD_E, or OSD_CS"))
        0 <= osd_order <= typemax(Int) ||
            throw(ArgumentError("osd_order must be a nonnegative Int"))
        return new(Int(max_iter), String(bp_method), Float64(ms_scaling_factor),
            String(schedule), Int(omp_thread_count), String(osd_method),
            Int(osd_order))
    end
end

_toric3d_settings_metadata(s::Toric3DDecoderSettings) = (
    max_iter=s.max_iter, bp_method=s.bp_method,
    ms_scaling_factor=s.ms_scaling_factor, schedule=s.schedule,
    omp_thread_count=s.omp_thread_count, osd_method=s.osd_method,
    osd_order=s.osd_order)

"""GF(2) row reduction and kernel words for the deliberately small X oracle."""
struct Toric3DExactSolver
    transform::BitMatrix
    pivots::Vector{Int}
    kernel_words::Vector{BitVector}
end

function _toric3d_exact_solver(matrix::BitMatrix)
    m, n = size(matrix)
    reduced = copy(matrix)
    transform = falses(m, m)
    for i in 1:m
        transform[i, i] = true
    end
    pivots = Int[]
    pivot_row = 1
    for column in 1:n
        pivot_row > m && break
        candidate = findfirst(row -> reduced[row, column], pivot_row:m)
        candidate === nothing && continue
        selected = pivot_row + candidate - 1
        if selected != pivot_row
            reduced[pivot_row, :], reduced[selected, :] =
                copy(reduced[selected, :]), copy(reduced[pivot_row, :])
            transform[pivot_row, :], transform[selected, :] =
                copy(transform[selected, :]), copy(transform[pivot_row, :])
        end
        for row in 1:m
            if row != pivot_row && reduced[row, column]
                @views reduced[row, :] .= xor.(reduced[row, :], reduced[pivot_row, :])
                @views transform[row, :] .= xor.(transform[row, :], transform[pivot_row, :])
            end
        end
        push!(pivots, column)
        pivot_row += 1
    end

    free_columns = setdiff(collect(1:n), pivots)
    length(free_columns) <= 9 ||
        throw(ArgumentError("exact reference requires at most 512 kernel words"))
    words = BitVector[falses(n)]
    for free_column in free_columns
        basis = falses(n)
        basis[free_column] = true
        for (row, pivot) in enumerate(pivots)
            basis[pivot] = reduced[row, free_column]
        end
        append!(words, [word .⊻ basis for word in words])
    end
    return Toric3DExactSolver(transform, pivots, words)
end

function _toric3d_exact_decode(solver::Toric3DExactSolver,
        syndrome::AbstractVector{Bool}, n::Int)
    transformed = BitVector(isodd(count(solver.transform[row, :] .& syndrome))
        for row in axes(solver.transform, 1))
    rank = length(solver.pivots)
    any(@view transformed[rank+1:end]) &&
        throw(ArgumentError("X syndrome is inconsistent with the code"))
    particular = falses(n)
    for (row, pivot) in enumerate(solver.pivots)
        particular[pivot] = transformed[row]
    end
    best = particular
    best_weight = n + 1
    for word in solver.kernel_words
        candidate = particular .⊻ word
        weight = count(candidate)
        if weight < best_weight
            best, best_weight = candidate, weight
            best_weight == 0 && break
        end
    end
    return best
end

"""A reusable full-qubit CSS recovery decoder for one 3D toric slab."""
struct Toric3DDecoder
    code::ToricCode3D
    p::Float64
    x_backend::Symbol
    settings::Toric3DDecoderSettings
    star_matrix::BitMatrix
    plaquette_matrix::BitMatrix
    z_matcher::PythonCall.Py
    x_bposd::Union{Nothing,PythonCall.Py}
    x_exact::Union{Nothing,Toric3DExactSolver}
    z_version::String
    x_version::String
end

function _toric3d_numpy_matrix(numpy, matrix::BitMatrix)
    return numpy.asarray(Matrix{UInt8}(matrix); dtype=numpy.uint8)
end

function _toric3d_numpy_vector(numpy, bits::AbstractVector{Bool})
    return numpy.asarray(Vector{UInt8}(bits); dtype=numpy.uint8)
end

function _toric3d_py_bits(numpy, prediction, n::Int)
    contiguous = numpy.ascontiguousarray(prediction; dtype=numpy.uint8).reshape(-1)
    bytes = PythonCall.pyconvert(Vector{UInt8}, contiguous.tobytes())
    length(bytes) == n ||
        throw(ErrorException("decoder returned $(length(bytes)) qubits; expected $n"))
    all(bit -> bit <= 1, bytes) ||
        throw(ErrorException("decoder returned a nonbinary recovery"))
    return BitVector(bytes .== 1)
end

function _toric3d_load_ldpc()
    try
        return PythonCall.pyimport("ldpc")
    catch error
        throw(ErrorException(
            "ldpc BP+OSD could not be imported. Run `import CondaPkg; CondaPkg.resolve()` " *
            "for this project. Original error: $(sprint(showerror, error))"))
    end
end

function _toric3d_module_version(py_module, package::String)
    try
        return PythonCall.pyconvert(String, py_module.__version__)
    catch
        metadata = PythonCall.pyimport("importlib.metadata")
        return PythonCall.pyconvert(String, metadata.version(package))
    end
end

"""
    build_toric3d_decoder(code; p, x_backend=:bposd, settings=Toric3DDecoderSettings())

Build PyMatching Z recovery with unit weights and either ldpc BP+OSD or the
small-code exact minimum-data-weight X reference. The exact oracle supports
only `(Lx,Ly,Lz) = (1,2,2)` or `(2,2,2)`.
"""
function build_toric3d_decoder(code::ToricCode3D; p::Real,
        x_backend::Symbol=:bposd,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings())
    isfinite(p) && 0 <= p <= 1 ||
        throw(ArgumentError("p must lie in [0, 1]"))
    x_backend in (:bposd, :exact) ||
        throw(ArgumentError("x_backend must be :bposd or :exact"))
    if x_backend === :exact
        (code.Lx in (1, 2) && code.Ly == 2 && code.Lz == 2) ||
            throw(ArgumentError("exact reference is available only for (1 or 2,2,2)"))
    end
    star_matrix = star_check_matrix(code)
    plaquette_matrix = plaquette_check_matrix(code)
    all(vec(sum(star_matrix; dims=1)) .<= 2) ||
        throw(ArgumentError("star check matrix is not graphlike"))

    numpy = PythonCall.pyimport("numpy")
    pymatching = _load_pymatching()
    n = data_qubit_count(code)
    # PyMatching merges equal parallel columns at Ly/Lz=2. Default identity
    # fault IDs preserve one real qubit representative per edge; padding makes
    # the returned data-recovery vector length n even if final IDs were merged.
    z_matcher = pymatching.Matching.from_check_matrix(
        _toric3d_numpy_matrix(numpy, star_matrix);
        weights=1.0, merge_strategy="smallest-weight", min_num_fault_ids=n)
    z_version = _toric3d_module_version(pymatching, "pymatching")

    if x_backend === :exact
        exact = _toric3d_exact_solver(plaquette_matrix)
        return Toric3DDecoder(code, Float64(p), x_backend, settings,
            star_matrix, plaquette_matrix, z_matcher, nothing, exact,
            z_version, "builtin")
    end

    ldpc = _toric3d_load_ldpc()
    # p=0 is common for noiseless scans, but a manual nonzero syndrome must
    # still be decodable. A tiny positive prior avoids infinite BP log odds.
    prior = clamp(Float64(p), 1e-12, 1 - 1e-12)
    bposd = ldpc.BpOsdDecoder(_toric3d_numpy_matrix(numpy, plaquette_matrix);
        error_rate=prior, max_iter=settings.max_iter,
        bp_method=settings.bp_method,
        ms_scaling_factor=settings.ms_scaling_factor,
        schedule=settings.schedule,
        omp_thread_count=settings.omp_thread_count,
        osd_method=settings.osd_method,
        osd_order=settings.osd_order)
    return Toric3DDecoder(code, Float64(p), x_backend, settings,
        star_matrix, plaquette_matrix, z_matcher, bposd, nothing,
        z_version, _toric3d_module_version(ldpc, "ldpc"))
end

function _toric3d_matrix_syndrome(matrix::BitMatrix, frame::BitVector)
    length(frame) == size(matrix, 2) ||
        throw(DimensionMismatch("syndrome matrix and frame dimensions differ"))
    result = falses(size(matrix, 1))
    for column in findall(frame)
        # A copied BitVector column retains packed-word XOR; a view does not.
        result .⊻= matrix[:, column]
    end
    return result
end

"""Return a full-qubit Pauli recovery, or fail if either syndrome disagrees."""
function decode_toric3d(decoder::Toric3DDecoder,
        stars::AbstractVector{Bool}, plaquettes::AbstractVector{Bool})
    length(stars) == size(decoder.star_matrix, 1) ||
        throw(DimensionMismatch("star syndrome has the wrong length"))
    length(plaquettes) == size(decoder.plaquette_matrix, 1) ||
        throw(DimensionMismatch("plaquette syndrome has the wrong length"))
    n = data_qubit_count(decoder.code)
    !any(stars) && !any(plaquettes) &&
        return PauliFrame(falses(n), falses(n))

    numpy = PythonCall.pyimport("numpy")
    z = if any(stars)
        _toric3d_py_bits(numpy,
            decoder.z_matcher.decode(_toric3d_numpy_vector(numpy, stars)), n)
    else
        falses(n)
    end
    x = if !any(plaquettes)
        falses(n)
    elseif decoder.x_backend === :exact
        _toric3d_exact_decode(decoder.x_exact::Toric3DExactSolver, plaquettes, n)
    else
        _toric3d_py_bits(numpy,
            (decoder.x_bposd::PythonCall.Py).decode(
                _toric3d_numpy_vector(numpy, plaquettes)), n)
    end
    _toric3d_matrix_syndrome(decoder.star_matrix, z) == stars ||
        throw(ErrorException("Z recovery does not reproduce star syndrome"))
    _toric3d_matrix_syndrome(decoder.plaquette_matrix, x) == plaquettes ||
        throw(ErrorException("X recovery does not reproduce plaquette syndrome"))
    return PauliFrame(x, z)
end

"""Backend and reproducibility information suitable for scan outputs."""
toric3d_decoder_metadata(decoder::Toric3DDecoder) = (
    x_backend=decoder.x_backend === :bposd ? "ldpc BP+OSD" : "exact GF(2)",
    x_version=decoder.x_version,
    z_backend="PyMatching",
    z_version=decoder.z_version,
    p=decoder.p,
    settings=_toric3d_settings_metadata(decoder.settings),
)
