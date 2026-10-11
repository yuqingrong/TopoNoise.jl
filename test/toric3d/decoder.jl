using Test

function _toric3d_test_syndrome(matrix, error)
    return BitVector(isodd(count((@view matrix[row, :]) .& error))
                     for row in axes(matrix, 1))
end

@testset "3D toric full recovery decoders" begin
    small = ToricCode3D(2, 2, 2)
    hx = star_check_matrix(small)
    hz = plaquette_check_matrix(small)
    n = data_qubit_count(small)
    @test n == 28
    @test maximum(vec(sum(hx; dims=1))) <= 2

    @testset "exact reference minimizes X data weight" begin
        decoder = build_toric3d_decoder(small; p=0.0, x_backend=:exact)
        @test decode_toric3d(decoder, falses(size(hx, 1)), falses(size(hz, 1))) ==
              PauliFrame(falses(n), falses(n))

        # Every syndrome made by up to two X faults has a representative of
        # weight at most two. The exact oracle must find that minimum, even
        # when another representative has a different logical parity.
        for indices in ([i] for i in 1:n)
            x = falses(n)
            x[indices] .= true
            plaquettes = _toric3d_test_syndrome(hz, x)
            recovery = decode_toric3d(decoder, falses(size(hx, 1)), plaquettes)
            @test recovery.z == falses(n)
            @test _toric3d_test_syndrome(hz, recovery.x) == plaquettes
            @test count(recovery.x) <= 1
        end

        for i in 1:n-1, j in i+1:n
            x = falses(n)
            x[[i, j]] .= true
            plaquettes = _toric3d_test_syndrome(hz, x)
            recovery = decode_toric3d(decoder, falses(size(hx, 1)), plaquettes)
            @test _toric3d_test_syndrome(hz, recovery.x) == plaquettes
            @test count(recovery.x) <= 2
        end

        impossible = falses(size(hz, 1))
        impossible[1] = true
        @test_throws ArgumentError decode_toric3d(decoder,
            falses(size(hx, 1)), impossible)
    end

    @testset "matching returns full Z-qubit recovery with parallel columns" begin
        decoder = build_toric3d_decoder(small; p=0.08, x_backend=:exact)
        for i in 1:n
            z = falses(n)
            z[i] = true
            stars = _toric3d_test_syndrome(hx, z)
            recovery = decode_toric3d(decoder, stars, falses(size(hz, 1)))
            @test length(recovery.z) == n
            @test _toric3d_test_syndrome(hx, recovery.z) == stars
            @test count(recovery.z) <= 1
        end
    end

    @testset "BP+OSD settings, syndrome, and exact comparator" begin
        settings = Toric3DDecoderSettings(max_iter=17, bp_method="minimum_sum",
            ms_scaling_factor=0.75, schedule="parallel", omp_thread_count=1,
            osd_method="OSD_CS", osd_order=2)
        bp = build_toric3d_decoder(small; p=0.06, settings)
        exact = build_toric3d_decoder(small; p=0.06, x_backend=:exact, settings)
        meta = toric3d_decoder_metadata(bp)
        @test meta isa NamedTuple
        @test meta.x_backend == "ldpc BP+OSD"
        @test meta.z_backend == "PyMatching"
        @test !isempty(meta.x_version)
        @test !isempty(meta.z_version)
        @test meta.settings.max_iter == 17
        @test meta.settings.ms_scaling_factor == 0.75
        @test meta.settings.osd_order == 2

        zero_prior = build_toric3d_decoder(small; p=0.0, settings)
        one_fault = falses(n)
        one_fault[1] = true
        zero_prior_result = decode_toric3d(zero_prior, falses(size(hx, 1)),
            _toric3d_test_syndrome(hz, one_fault))
        @test _toric3d_test_syndrome(hz, zero_prior_result.x) ==
              _toric3d_test_syndrome(hz, one_fault)

        for indices in ([1], [2, 7], [4, 12, 19], [3, 10, 25])
            x = falses(n)
            x[indices] .= true
            plaquettes = _toric3d_test_syndrome(hz, x)
            bp_recovery = decode_toric3d(bp, falses(size(hx, 1)), plaquettes)
            exact_recovery = decode_toric3d(exact, falses(size(hx, 1)), plaquettes)
            @test _toric3d_test_syndrome(hz, bp_recovery.x) == plaquettes
            @test count(exact_recovery.x) <= count(bp_recovery.x)
        end

        @test_throws DimensionMismatch decode_toric3d(bp, falses(size(hx, 1) - 1),
            falses(size(hz, 1)))
        @test_throws DimensionMismatch decode_toric3d(bp, falses(size(hx, 1)),
            falses(size(hz, 1) - 1))
        @test_throws ArgumentError build_toric3d_decoder(small; p=-0.1)
        @test_throws ArgumentError build_toric3d_decoder(small; p=1.1)
        @test_throws ArgumentError build_toric3d_decoder(small; p=0.1, x_backend=:bogus)
        @test_throws ArgumentError build_toric3d_decoder(ToricCode3D(3, 2, 2);
            p=0.1, x_backend=:exact)
    end
end
