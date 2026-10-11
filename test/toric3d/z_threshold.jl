using Test, Random

@testset "Logical Z threshold API" begin
    @test isdefined(TopoNoise, :build_toric3d_z_fault_model)
    @test isdefined(TopoNoise, :decode_toric3d_z_batch)
    @test isdefined(TopoNoise, :scan_toric3d_z_threshold)
end

if isdefined(TopoNoise, :build_toric3d_z_fault_model)
    @testset "Exact Z fault map agrees with every explicit gate fault" begin
        code = Toric3DEncodedSlice(2,2,3)
        model = build_toric3d_z_fault_model(code)
        encoder = model.encoder
        clean = sample_toric3d_faults(MersenneTwister(2),encoder)
        @test size(model.matrix) == (data_qubit_count(code),3*2*3+13*2*2*3)
        @test all(1 .<= vec(sum(model.matrix;dims=1)) .<= 3)
        @test all(iszero, rem.(sum(star_check_matrix(code)*model.matrix;dims=1),2))
        for (column,location) in enumerate(model.locations)
            expected = falses(data_qubit_count(code))
            if location.kind === :seed
                expected[location.qubit] = true
            else
                record = with_toric3d_fault(encoder,clean,location.step,location.qubit,:Z)
                frame = PauliFrame(copy(expected),copy(expected))
                propagate_toric3d_segment!(frame,encoder,record,1,length(encoder.operations))
                expected .= frame.z
            end
            actual = BitVector(model.matrix[:,column] .!= 0)
            @test actual == expected
            @test count(isodd.(star_check_matrix(code)*actual)) == 2
        end
    end

    @testset "Encoded MWPM validates full data recovery" begin
        code = Toric3DEncodedSlice(2,5,5)
        model = build_toric3d_z_fault_model(code)
        decoder = build_toric3d_z_threshold_decoder(code)
        frames = BitMatrix(Matrix(model.matrix) .!= 0)
        counts = toric3d_z_threshold_counts(decoder,frames)
        @test counts.failures == counts.logical_1_failures == counts.logical_2_failures == 0
        H = star_check_matrix(code)
        syndromes = BitMatrix(isodd.(H*frames))
        recovery = decode_toric3d_z_batch(decoder,syndromes)
        @test isodd.(H*recovery) == syndromes
        bad = copy(recovery); bad[1,1] ⊻= true
        @test_throws ErrorException TopoNoise._toric3d_validate_z_batch(decoder,syndromes,bad)
        @test_throws DimensionMismatch decode_toric3d_z_batch(decoder,falses(2,3))
        odd = falses(size(H,1),1); odd[1] = true
        @test_throws ArgumentError decode_toric3d_z_batch(decoder,odd)
        for (index,support) in enumerate(toric3d_encoded_logical_z(code))
            logical = falses(data_qubit_count(code),1); logical[support,1] .= true
            count = toric3d_z_threshold_counts(decoder,logical)
            @test count.failures == 1
            @test getproperty(count,Symbol("logical_",index,"_failures")) == 1
        end
        trivial = falses(data_qubit_count(code),1)
        trivial[first(plaquette_checks(code)),1] .= true
        @test toric3d_z_threshold_counts(decoder,trivial).failures == 0
    end

    @testset "Z threshold sampling and incremental scan" begin
        model = build_toric3d_z_fault_model(Toric3DEncodedSlice(1,2,3))
        a = sample_toric3d_z_frames(MersenneTwister(17),model;p=.03,shots=12)
        b = sample_toric3d_z_frames(MersenneTwister(17),model;p=.03,shots=12)
        @test a == b
        @test !any(sample_toric3d_z_frames(MersenneTwister(17),model;p=0,shots=4))
        @test_throws ArgumentError sample_toric3d_z_frames(MersenneTwister(1),model;p=-.1,shots=2)
        @test_throws ArgumentError sample_toric3d_z_frames(MersenneTwister(1),model;p=.1,shots=0)
        args = (sizes=[(2,2,3)],error_rates=[0.,.01],shots=20,seed=19,batch_size=7)
        mktempdir() do directory
            a = scan_toric3d_z_threshold(MersenneTwister(19);args...,output_dir=directory)
            b = scan_toric3d_z_threshold(MersenneTwister(19);args...)
            @test length(a.points) == 2
            @test [p.failures for p in a.points] == [p.failures for p in b.points]
            @test first(a.points).failures == 0
            @test isfile(joinpath(directory,"toric3d-z-threshold-raw.csv"))
            @test length(readlines(joinpath(directory,"toric3d-z-threshold-raw.csv"))) == 3
            metadata = TopoNoise.TOML.parsefile(joinpath(directory,"toric3d-z-threshold-metadata.toml"))
            @test metadata["complete"]
            @test metadata["seed"] == 19
            @test metadata["logical_qubits"] == 2
        end
    end
end
