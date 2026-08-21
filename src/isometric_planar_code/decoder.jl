function _isometric_planar_matching(model::IsometricPlanarCodeModel)
    distance = model.distance
    distance == 2 && throw(ArgumentError(
        "distance-2 planar decoding uses its deterministic diagnostic lookup"))
    return get!(_isometric_planar_matching_cache, distance) do
        matching = _pymatching().Matching()
        plaquette(row, col) =
            _isometric_planar_plaquette_detector(distance, row, col)
        north(col) = _isometric_planar_north_detector(distance, col)
        south(col) = _isometric_planar_south_detector(distance, col)
        logical_id = _isometric_planar_fault_count(distance)

        for row in 1:distance, col in 1:(distance - 1)
            fault = _isometric_planar_h_fault(distance, row, col)
            first, second = if row == 1
                north(col), plaquette(1, col)
            elseif row == distance
                plaquette(distance - 1, col), south(col)
            else
                plaquette(row - 1, col), plaquette(row, col)
            end
            matching.add_edge(pyint(first), pyint(second);
                fault_ids=pyset([fault]), weight=1.0,
                merge_strategy="independent")
        end

        for row in 1:(distance - 1), col in 1:distance
            fault = _isometric_planar_v_fault(distance, row, col)
            fault_ids = col == 1 ? [fault, logical_id] : [fault]
            if col == 1 || col == distance
                adjacent = col == 1 ? plaquette(row, 1) :
                                         plaquette(row, distance - 1)
                matching.add_boundary_edge(pyint(adjacent);
                    fault_ids=pyset(fault_ids), weight=1.0,
                    merge_strategy="independent")
            else
                matching.add_edge(pyint(plaquette(row, col - 1)),
                    pyint(plaquette(row, col)); fault_ids=pyset(fault_ids),
                    weight=1.0, merge_strategy="independent")
            end
        end
        matching
    end
end

function _validate_isometric_planar_syndrome(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    expected_plaquettes = (model.distance - 1, model.distance - 1)
    size(syndrome.plaquettes) == expected_plaquettes || throw(DimensionMismatch(
        "plaquette syndrome must have size $expected_plaquettes"))
    expected_boundary = model.distance - 1
    length(syndrome.north) == expected_boundary || throw(DimensionMismatch(
        "north syndrome must have length $expected_boundary"))
    length(syndrome.south) == expected_boundary || throw(DimensionMismatch(
        "south syndrome must have length $expected_boundary"))
    return nothing
end

function _isometric_planar_syndrome_vector(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    _validate_isometric_planar_syndrome(model, syndrome)
    distance = model.distance
    values = zeros(Int8, _isometric_planar_detector_count(distance))
    for row in 1:(distance - 1), col in 1:(distance - 1)
        detector = _isometric_planar_plaquette_detector(distance, row, col)
        values[detector + 1] = syndrome.plaquettes[row, col]
    end
    for col in 1:(distance - 1)
        values[_isometric_planar_north_detector(distance, col) + 1] = syndrome.north[col]
        values[_isometric_planar_south_detector(distance, col) + 1] = syndrome.south[col]
    end
    return values
end

"""
    decode_isometric_planar_syndrome(model, syndrome) -> Correction

Decode plaquette and local north/south detector bits. Distance two uses a
fixed west-boundary representative for its unavoidable one-fault tie and is
provided only as a diagnostic geometry.
"""
function decode_isometric_planar_syndrome(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    syndrome_vector = _isometric_planar_syndrome_vector(model, syndrome)
    correction = Correction(model.rows, model.cols)
    any(value -> !iszero(value), syndrome_vector) || return correction
    if model.distance == 2
        correction.vertical[1, 1] = true
        return correction
    end

    prediction = _isometric_planar_matching(model).decode(
        pyimport("numpy").asarray(syndrome_vector))
    values = pyconvert(Vector{Int}, prediction)
    fault_count = _isometric_planar_fault_count(model.distance)
    length(values) >= fault_count || throw(ErrorException(
        "PyMatching returned $(length(values)) fault predictions; expected at least $fault_count"))
    for row in 1:model.rows, col in 1:(model.cols - 1)
        correction.horizontal[row, col] =
            values[_isometric_planar_h_fault(model.distance, row, col) + 1] == 1
    end
    for row in 1:(model.rows - 1), col in 1:model.cols
        correction.vertical[row, col] =
            values[_isometric_planar_v_fault(model.distance, row, col) + 1] == 1
    end
    return correction
end

"""Return whether the error/correction residual has odd west-boundary parity."""
function isometric_planar_logical_failure(
        errors::IsometricPlanarVirtualErrors, correction::Correction)
    size(errors.horizontal) == size(correction.horizontal) || throw(DimensionMismatch(
        "correction.horizontal must match horizontal virtual links"))
    size(errors.vertical) == size(correction.vertical) || throw(DimensionMismatch(
        "correction.vertical must match vertical virtual links"))
    parity = false
    for row in axes(errors.vertical, 1)
        parity = xor(parity, errors.vertical[row, 1], correction.vertical[row, 1])
    end
    return parity
end

function _isometric_planar_minimum_odd_logical_weight(model::IsometricPlanarCodeModel)
    distance = model.distance
    detector_count = _isometric_planar_detector_count(distance)
    boundary = detector_count + 1
    adjacency = [Tuple{Int,Bool}[] for _ in 1:boundary]
    add_edge!(first::Int, second::Int, logical::Bool=false) = begin
        push!(adjacency[first], (second, logical))
        push!(adjacency[second], (first, logical))
    end
    plaquette(row, col) =
        _isometric_planar_plaquette_detector(distance, row, col) + 1
    north(col) = _isometric_planar_north_detector(distance, col) + 1
    south(col) = _isometric_planar_south_detector(distance, col) + 1

    for row in 1:distance, col in 1:(distance - 1)
        row == 1 ? add_edge!(north(col), plaquette(1, col)) :
        row == distance ? add_edge!(plaquette(distance - 1, col), south(col)) :
        add_edge!(plaquette(row - 1, col), plaquette(row, col))
    end
    for row in 1:(distance - 1), col in 1:distance
        col == 1 ? add_edge!(plaquette(row, 1), boundary, true) :
        col == distance ? add_edge!(plaquette(row, distance - 1), boundary) :
        add_edge!(plaquette(row, col - 1), plaquette(row, col))
    end

    distances = fill(-1, boundary, 2)
    queue = Tuple{Int,Bool}[(boundary, false)]
    distances[boundary, 1] = 0
    index = 1
    while index <= length(queue)
        node, parity = queue[index]
        index += 1
        node == boundary && parity && return distances[node, 2]
        for (next_node, edge_parity) in adjacency[node]
            next_parity = xor(parity, edge_parity)
            parity_index = next_parity ? 2 : 1
            distances[next_node, parity_index] >= 0 && continue
            distances[next_node, parity_index] = distances[node, parity ? 2 : 1] + 1
            push!(queue, (next_node, next_parity))
        end
    end
    throw(ErrorException("planar detector graph has no odd-logical path"))
end

function _isometric_planar_detector_batch(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel,
        error_rate::Float64, shots::Int)
    distance = model.distance
    horizontal = rand(rng, distance, distance - 1, shots) .< error_rate
    vertical = rand(rng, distance - 1, distance, shots) .< error_rate
    detectors = falses(shots, _isometric_planar_detector_count(distance))
    actual = falses(shots)
    for shot in 1:shots
        for row in 1:(distance - 1), col in 1:(distance - 1)
            detector = _isometric_planar_plaquette_detector(distance, row, col)
            detectors[shot, detector + 1] = xor(
                horizontal[row, col, shot], horizontal[row + 1, col, shot],
                vertical[row, col, shot], vertical[row, col + 1, shot])
        end
        for col in 1:(distance - 1)
            detectors[shot, _isometric_planar_north_detector(distance, col) + 1] =
                horizontal[1, col, shot]
            detectors[shot, _isometric_planar_south_detector(distance, col) + 1] =
                horizontal[distance, col, shot]
        end
        parity = false
        for row in 1:(distance - 1)
            parity = xor(parity, vertical[row, 1, shot])
        end
        actual[shot] = parity
    end
    return detectors, actual
end

function _isometric_planar_batch_predictions(
        model::IsometricPlanarCodeModel, detectors::BitMatrix)
    shots = size(detectors, 1)
    size(detectors, 2) == _isometric_planar_detector_count(model.distance) ||
        throw(DimensionMismatch("detector batch has the wrong detector count"))
    if model.distance == 2
        return BitVector(any(detectors[shot, :]) for shot in 1:shots)
    end
    prediction = _isometric_planar_matching(model).decode_batch(
        pyimport("numpy").asarray(detectors))
    values = pyconvert(BitMatrix, prediction)
    expected = _isometric_planar_fault_count(model.distance) + 1
    size(values) == (shots, expected) || throw(ErrorException(
        "PyMatching returned predictions with size $(size(values)); expected ($shots, $expected)"))
    return BitVector(values[:, end])
end

"""One batched logical-failure estimate for the planar virtual-bond model."""
