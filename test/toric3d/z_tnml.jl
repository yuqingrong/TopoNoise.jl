using Test
using TopoNoise
using SparseArrays
using TOML

@testset "Z-only TN package API" begin
    for name in (:Toric3DZTNSettings, :Toric3DZTNTopology, :Toric3DZTNDecoder,
            :build_toric3d_z_tn_topology, :toric3d_z_tn_edge_probabilities,
            :build_toric3d_z_tn_decoder, :decode_toric3d_z_tn)
        @test isdefined(TopoNoise, name)
        @test name in names(TopoNoise)
    end
end

# Enumerate independent faults directly. This oracle neither contracts tensors
# nor uses the decoder's edge-grouping or normalization helpers.
function z_tnml_enumerated_weights(n, edges, probabilities, logical_masks)
    weights = zeros(4, 1 << n)
    for fault_mask in 0:((1 << length(edges)) - 1)
        syndrome_mask = 0
        logical_mask = 0
        probability = 1.0
        for (j, ((u, v), q, logical)) in enumerate(zip(edges, probabilities, logical_masks))
            active = !iszero(fault_mask & (1 << (j - 1)))
            probability *= active ? q : 1 - q
            if active
                syndrome_mask ⊻= (1 << (u - 1)) ⊻ (1 << (v - 1))
                logical_mask ⊻= Int(logical)
            end
        end
        weights[logical_mask + 1, syndrome_mask + 1] += probability
    end
    return weights
end

function z_tnml_syndrome(n, mask)
    return BitVector(!iszero(mask & (1 << (j - 1))) for j in 1:n)
end

@testset "Validated defaults and rejected numerical controls" begin
    settings = Toric3DZTNSettings()
    @test settings.max_bond_dimension == 240
    @test settings.cutoff == 1e-18
    @test settings.negative_tolerance == 1e-10
    for dimension in (0, -1)
        @test_throws ArgumentError Toric3DZTNSettings(max_bond_dimension=dimension)
    end
    for bad in (-1.0, 1.0, Inf, NaN)
        @test_throws ArgumentError Toric3DZTNSettings(cutoff=bad)
        @test_throws ArgumentError Toric3DZTNSettings(negative_tolerance=bad)
    end
end

@testset "Reduced topology agrees with individual propagated circuit faults" begin
    code = Toric3DEncodedSlice(1, 2, 2)
    model = build_toric3d_z_fault_model(code)
    topology = build_toric3d_z_tn_topology(model)
    checks = star_checks(code)
    logicals = toric3d_encoded_logical_x(code)

    # Compute each signature from the final physical support, independently of
    # the matrix multiplication and signature grouping in the topology builder.
    for column in axes(model.matrix, 2)
        support = BitVector(!iszero(x) for x in model.matrix[:, column])
        expected_detectors = BitVector(isodd(count(support[check])) for check in checks)
        expected_logicals = BitVector(isodd(count(support[check])) for check in logicals)
        @test topology.raw_detector_matrix[:, column] == UInt8.(expected_detectors)
        @test topology.raw_logical_matrix[:, column] == UInt8.(expected_logicals)
    end
    @test all(==(2), vec(sum(topology.raw_detector_matrix; dims=1)))
    @test allunique(collect(zip(topology.detector_endpoints, topology.logical_masks)))
    @test topology.multiplicities == length.(topology.original_columns)
    @test sort!(reduce(vcat, topology.original_columns)) == collect(axes(model.matrix, 2))
    for (edge, mask, columns) in zip(topology.detector_endpoints,
            topology.logical_masks, topology.original_columns), column in columns
        @test Tuple(findall(!iszero, topology.raw_detector_matrix[:, column])) == edge
        @test topology.raw_logical_matrix[1, column] +
              2topology.raw_logical_matrix[2, column] == mask
    end
end

@testset "Parallel fault reduction preserves detector and logical signatures" begin
    # Columns 1 and 2 may be reduced together. Column 3 has the same endpoints
    # but a different logical bit, so it must remain a separate effective edge.
    detectors = sparse(UInt8[1 1 1 0 1; 1 1 1 1 0; 0 0 0 1 1])
    logicals = sparse(UInt8[0 0 1 0 1; 0 0 0 1 1])
    topology = Toric3DZTNTopology(detectors, logicals,
        [(1, 2), (1, 2), (2, 3), (1, 3)], UInt8[0, 1, 2, 3],
        [2, 1, 1, 1], [[1, 2], [3], [4], [5]])
    p = 0.13
    @test toric3d_z_tn_edge_probabilities(topology, p) ≈
          [2p * (1 - p), p, p, p] atol=1e-15
    @test toric3d_z_tn_edge_probabilities(topology, 0.0) == zeros(4)
    @test toric3d_z_tn_edge_probabilities(topology, 1e-20) ≈
          [2e-20, 1e-20, 1e-20, 1e-20] rtol=1e-14

    raw_edges = [(1, 2), (1, 2), (1, 2), (2, 3), (1, 3)]
    exact = z_tnml_enumerated_weights(3, raw_edges, fill(p, 5), UInt8[0, 0, 1, 2, 3])
    decoder = build_toric3d_z_tn_decoder(topology, p;
        settings=Toric3DZTNSettings(max_bond_dimension=128, cutoff=0.0))
    @test decoder isa Toric3DZTNDecoder
    for mask in 0:7
        isodd(count_ones(mask)) && continue
        expected = exact[:, mask + 1]
        result = decode_toric3d_z_tn(decoder, z_tnml_syndrome(3, mask))
        @test result.valid
        @test result.probabilities ≈ expected / sum(expected) atol=2e-12 rtol=2e-12
        @test result.log_syndrome ≈ log(sum(expected)) atol=2e-12 rtol=2e-12
        @test result.logical_mask == argmax(expected) - 1
    end

    for p_bad in (-0.01, 0.5, 1.0, Inf, NaN)
        @test_throws ArgumentError toric3d_z_tn_edge_probabilities(topology, p_bad)
        @test_throws ArgumentError build_toric3d_z_tn_decoder(topology, p_bad)
    end
    @test_throws DimensionMismatch Toric3DZTNTopology(detectors, logicals,
        [(1, 2)], UInt8[0, 1, 2, 3], [2, 1, 1, 1], [[1, 2], [3], [4], [5]])
    @test_throws ArgumentError Toric3DZTNTopology(detectors, logicals,
        [(1, 2), (1, 2), (2, 3), (1, 3)], UInt8[0, 0, 2, 3],
        [2, 1, 1, 1], [[1, 2], [3], [4], [5]])
    @test_throws ArgumentError Toric3DZTNTopology(detectors, logicals,
        [(1, 2), (1, 2), (2, 3), (1, 3)], UInt8[0, 1, 2, 3],
        [2, 1, 1, 1], [[1, 2], [2], [4], [5]])
end

@testset "Conditioned-star contraction agrees with independent enumeration" begin
    fixtures = [
        (4, [(1, 2), (2, 3), (3, 4), (1, 4), (1, 2)],
            [0.03, 0.11, 0.07, 0.13, 0.05], UInt8[0, 1, 2, 3, 1]),
        (5, [(1, 2), (1, 3), (1, 4), (1, 5)],
            [0.1, 0.2, 0.12, 0.08], UInt8[0, 1, 2, 3]),
        (4, [(1, 2), (2, 3), (3, 4)],
            [0.04, 0.09, 0.15], UInt8[1, 2, 3]),
    ]
    settings = Toric3DZTNSettings(max_bond_dimension=128, cutoff=0.0)
    for (n, edges, qs, logicals) in fixtures
        exact = z_tnml_enumerated_weights(n, edges, qs, logicals)
        decoder = TopoNoise._build_toric3d_z_tn_edges(n, edges, qs, logicals; settings)
        for mask in 0:((1 << n) - 1)
            isodd(count_ones(mask)) && continue
            expected = exact[:, mask + 1]
            result = decode_toric3d_z_tn(decoder, z_tnml_syndrome(n, mask))
            @test result.valid
            @test result.probabilities ≈ expected / sum(expected) atol=2e-12 rtol=2e-12
            @test result.log_syndrome ≈ log(sum(expected)) atol=2e-12 rtol=2e-12
            @test result.logical_mask == argmax(expected) - 1
            @test result.raw_probabilities == result.probabilities
            @test result.diagnostics.roundoff_correction_l1 == 0.0
            @test result.diagnostics.maximum_local_truncation_error <= 1e-25
        end
    end
end

@testset "Conditioned-star truncation visits only the active prefix" begin
    # The first graph is a chain. The second has gaps and two untouched tail
    # vertices, so it catches implementations that sweep future detector sites.
    fixtures = [
        (5, [(1, 2), (2, 3), (3, 4), (4, 5)],
            [0.03, 0.07, 0.11, 0.13], UInt8[0, 1, 2, 3], 4),
        (8, [(1, 4), (1, 6), (3, 5), (5, 6)],
            [0.04, 0.09, 0.12, 0.08], UInt8[1, 2, 3, 0], 15),
    ]
    settings = Toric3DZTNSettings(max_bond_dimension=128, cutoff=0.0)
    for (n, edges, qs, logicals, expected_events) in fixtures
        exact = z_tnml_enumerated_weights(n, edges, qs, logicals)
        decoder = TopoNoise._build_toric3d_z_tn_edges(n, edges, qs, logicals; settings)
        for mask in 0:((1 << n) - 1)
            isodd(count_ones(mask)) && continue
            expected = exact[:, mask + 1]
            iszero(sum(expected)) && continue
            result = decode_toric3d_z_tn(decoder, z_tnml_syndrome(n, mask))
            @test result.valid
            @test result.probabilities ≈ expected / sum(expected) atol=2e-12 rtol=2e-12
            @test result.log_syndrome ≈ log(sum(expected)) atol=2e-12 rtol=2e-12
            @test result.diagnostics.events == expected_events
            @test result.diagnostics.local_svd_events == expected_events
        end
    end
end

@testset "Logical-mask bit order and parallel-edge parity" begin
    settings = Toric3DZTNSettings(cutoff=0.0)
    for logical_mask in 0:3
        decoder = TopoNoise._build_toric3d_z_tn_edges(2, [(1, 2)], [0.2],
            UInt8[logical_mask]; settings)
        result = decode_toric3d_z_tn(decoder, trues(2))
        @test result.probabilities ≈ [i == logical_mask + 1 ? 1.0 : 0.0 for i in 1:4] atol=1e-14
        @test result.logical_mask == logical_mask
    end
    decoder = TopoNoise._build_toric3d_z_tn_edges(2, [(1, 2), (1, 2)],
        [0.5, 0.5], UInt8[1, 2]; settings)
    result = decode_toric3d_z_tn(decoder, falses(2))
    @test result.probabilities ≈ [0.5, 0.0, 0.0, 0.5] atol=1e-14
    @test result.logical_mask in (0, 3)
end

@testset "Zero noise and invalid decoder inputs through the public API" begin
    code = Toric3DEncodedSlice(1, 2, 2)
    model = build_toric3d_z_fault_model(code)
    topology = build_toric3d_z_tn_topology(model)
    n = size(topology.raw_detector_matrix, 1)
    for source in (code, model, topology)
        decoder = build_toric3d_z_tn_decoder(source, 0.0)
        result = decode_toric3d_z_tn(decoder, falses(n))
        @test result.valid
        @test result.probabilities == [1.0, 0.0, 0.0, 0.0]
        @test result.raw_probabilities == result.probabilities
        @test result.logical_mask == 0
        @test result.log_syndrome == 0.0
        @test result.diagnostics.detector_count == n
        @test result.diagnostics.edge_count == length(topology.detector_endpoints)
        @test_throws DimensionMismatch decode_toric3d_z_tn(decoder, falses(n - 1))
        odd = falses(n); odd[1] = true
        @test_throws ArgumentError decode_toric3d_z_tn(decoder, odd)
        nonbinary = zeros(Int, n); nonbinary[1] = 2
        @test_throws ArgumentError decode_toric3d_z_tn(decoder, nonbinary)
        impossible = falses(n); impossible[1:2] .= true
        @test_throws DomainError decode_toric3d_z_tn(decoder, impossible)
    end
    for source in (code, model), p_bad in (-0.01, 0.5, Inf, NaN)
        @test_throws ArgumentError build_toric3d_z_tn_decoder(source, p_bad)
    end
end

@testset "Explicit truncation reports local diagnostics and signed weights" begin
    edges = [(u, v) for u in 1:5 for v in (u + 1):5]
    settings = Toric3DZTNSettings(max_bond_dimension=2, cutoff=0.0)
    decoder = TopoNoise._build_toric3d_z_tn_edges(5, edges, fill(0.1, length(edges)),
        UInt8.(mod.(1:length(edges), 4)); settings)
    result = decode_toric3d_z_tn(decoder, falses(5))
    @test result.diagnostics.actual_maxbond <= 2
    @test result.diagnostics.events > 0
    @test result.diagnostics.discarded_weight > 0
    @test all(isfinite, result.probabilities)
    @test sum(result.probabilities) ≈ 1.0 atol=1e-14
    @test result.probabilities == result.raw_probabilities
    @test result.diagnostics.roundoff_correction_l1 == 0.0
    @test result.valid == (minimum(result.probabilities) >= -settings.negative_tolerance)
    @test result.logical_mask == (result.valid ? argmax(result.probabilities) - 1 : -1)
    # These are local SVD diagnostics; no posterior-error bound is inferred.
    @test result.diagnostics.sum_local_truncation_error >=
          result.diagnostics.maximum_local_truncation_error >= 0.0
end

@testset "Frozen L3 public-API regression against validated and exact posteriors" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "fixtures", "z_tnml_l3.toml"))
    code = Toric3DEncodedSlice(fixture["L"], fixture["L"], fixture["L"])
    model = build_toric3d_z_fault_model(code)
    topology = build_toric3d_z_tn_topology(model)
    @test size(topology.raw_detector_matrix) == (45, 378)
    @test length(topology.detector_endpoints) == 162
    # Use the public code constructor and validated defaults for the expensive
    # regression; a fresh Julia process needs only `using TopoNoise`.
    decoder = build_toric3d_z_tn_decoder(code, 0.02)
    for case in fixture["cases"]
        @testset "$(case["id"])" begin
            syndrome = falses(45); syndrome[case["syndrome_ones"]] .= true
            raw_faults = zeros(UInt8, 378); raw_faults[case["raw_fault_ones"]] .= 1
            @test isodd.(topology.raw_detector_matrix * raw_faults) == syndrome
            logical_bits = isodd.(topology.raw_logical_matrix * raw_faults)
            @test Int(logical_bits[1]) + 2Int(logical_bits[2]) == case["true_mask"]
            result = decode_toric3d_z_tn(decoder, syndrome)
            @test result.valid
            @test result.probabilities ≈ case["validated_probabilities"] atol=1e-10 rtol=1e-10
            @test result.log_syndrome ≈ case["validated_log_syndrome"] atol=1e-9 rtol=0.0
            @test result.logical_mask == case["validated_logical_mask"] == case["exact_mask"]
            @test sum(abs.(result.probabilities .- case["exact_probabilities"])) <= 0.001
            @test abs(result.log_syndrome - case["exact_log_marginal"]) <= 0.001
            @test result.raw_probabilities == result.probabilities
            @test result.diagnostics.actual_maxbond <= 240
            @test result.diagnostics.local_svd_events == 714
            @test result.diagnostics.algorithm ==
                  "conditioned_star_active_boundary_mps"
            @test result.diagnostics.roundoff_correction_l1 == 0.0
        end
    end
end
