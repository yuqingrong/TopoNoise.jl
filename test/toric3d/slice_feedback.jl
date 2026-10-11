using Test, Random

@testset "Slice-feedback public API" begin
    @test isdefined(TopoNoise, :build_toric3d_slice_decoder)
    @test isdefined(TopoNoise, :decode_toric3d_slice)
    @test isdefined(TopoNoise, :replay_toric3d_slice)
end

if isdefined(TopoNoise, :replay_toric3d_slice)
    slice_parity(H, bits) = BitVector(isodd.(Int.(H) * Int.(bits)))

    @testset "Slice matching returns physical recoveries and local routes" begin
        for dims in ((2,2,2),(2,2,3),(2,4,5))
            code = Toric3DSliceLattice(dims...)
            decoder = build_toric3d_slice_decoder(code)
            Hq = toric3d_slice_check_matrix(code)
            for q in axes(Hq,2)
                syndrome = Hq[:,q]
                result = decode_toric3d_slice(decoder,syndrome)
                @test slice_parity(Hq,result.recovery) == syndrome
                @test slice_parity(Hq,result.route) == syndrome
                @test !any(result.residual_syndrome)
                @test count(result.route) == 1
            end
            @test_throws DimensionMismatch decode_toric3d_slice(decoder,falses(1))
            odd = falses(size(Hq,1)); odd[1] = true
            @test_throws ArgumentError decode_toric3d_slice(decoder,odd)
            @test_throws ArgumentError decode_toric3d_slice(decoder,falses(size(Hq,1));routing=:invalid)
        end
        code = Toric3DSliceLattice(2,9,9)
        decoder = build_toric3d_slice_decoder(code)
        syndrome = falses(81); syndrome[1] = syndrome[4*9+1] = true
        step = decode_toric3d_slice(decoder,syndrome)
        full = decode_toric3d_slice(decoder,syndrome;routing=:full_path)
        @test length(only(step.paths)) == 4
        @test count(step.route) == 2
        @test count(step.residual_syndrome) == 2
        @test !any(full.residual_syndrome)
        @test count(full.route) == 4
        @test length(only(decode_toric3d_slice(decoder,step.residual_syndrome).paths)) == 2
    end

    @testset "Target feed-forward equals virtual control repair" begin
        code = Toric3DSliceLattice(2,3,3)
        encoder = toric3d_slice_encoder(code)
        decoder = build_toric3d_slice_decoder(code)
        Hq = toric3d_slice_check_matrix(code)
        zero_faults = sample_toric3d_faults(MersenneTwister(1),encoder)
        n = data_qubit_count(code)
        old = toric3d_slice_qubits(code,0)
        new = toric3d_slice_qubits(code,1)
        for routing in (:one_step,:full_path), seed in 1:12
            b = bitrand(MersenneTwister(seed),length(old))
            route = decode_toric3d_slice(decoder,slice_parity(Hq,b);routing).route
            physical = PauliFrame(falses(n),falses(n)); physical.x[old] .= b
            virtual = PauliFrame(copy(physical.x),copy(physical.z))
            virtual.x[old] .⊻= route
            propagate_toric3d_segment!(physical,encoder,zero_faults,1,encoder.checkpoints[1])
            propagate_toric3d_segment!(virtual,encoder,zero_faults,1,encoder.checkpoints[1])
            delta = TopoNoise._toric3d_slice_feedforward(code,route)
            physical.x[new] .⊻= delta
            virtual.x[old] .⊻= route
            @test physical == virtual
            residual = slice_parity(Hq,b .⊻ route)
            expected = BitVector(residual[mod(u-1,code.Ly)*code.Lz+v+1]
                                 for u in 0:code.Ly-1 for v in 0:code.Lz-1)
            @test slice_parity(Hq,physical.x[new]) == expected
            @test physical.x[old] == b
        end
    end

    @testset "Noiseless slice preparation and seeded fault replay" begin
        encoder = toric3d_slice_encoder(Toric3DSliceLattice(3,2,3))
        decoder = build_toric3d_slice_decoder(encoder.code)
        clean = sample_toric3d_faults(MersenneTwister(1),encoder)
        for feedback in (false,true), routing in (:one_step,:full_path)
            result = replay_toric3d_slice(encoder,clean,decoder;feedback,routing,record_history=true)
            @test !any(result.frame.x) && !any(result.frame.z)
            @test all(!any(s) for s in values(result.syndromes))
            @test length(result.history) == encoder.code.steps
        end
        a = sample_toric3d_faults(MersenneTwister(1234),encoder;p_x=.03,p_z=.02)
        b = sample_toric3d_faults(MersenneTwister(1234),encoder;p_x=.03,p_z=.02)
        @test a == b
        @test replay_toric3d_slice(encoder,a,decoder).frame == replay_toric3d_slice(encoder,b,decoder).frame
        @test_throws ArgumentError replay_toric3d_slice(encoder,a,decoder;routing=:bad)
        other = build_toric3d_slice_decoder(Toric3DSliceLattice(2,2,3))
        @test_throws DimensionMismatch replay_toric3d_slice(encoder,a,other)
        @test_throws ArgumentError with_toric3d_fault(encoder,a,1,encoder.operations[2].qubits[2],:Z)
    end

    @testset "Every single Z output fault stays within two slices and weight three" begin
        encoder = toric3d_slice_encoder(Toric3DSliceLattice(5,2,2))
        decoder = build_toric3d_slice_decoder(encoder.code)
        clean = sample_toric3d_faults(MersenneTwister(1),encoder)
        block = 3encoder.code.Ly*encoder.code.Lz
        for (gate,op) in enumerate(encoder.operations), q in op.qubits
            record = with_toric3d_fault(encoder,clean,gate,q,:Z)
            a = replay_toric3d_slice(encoder,record,decoder;feedback=false)
            b = replay_toric3d_slice(encoder,record,decoder;feedback=true)
            support = findall(a.frame.z)
            times = (support .- 1) .÷ block
            @test 1 <= length(support) <= 3
            @test maximum(times)-minimum(times) <= 1
            @test a.frame == b.frame
            @test !any(a.syndromes.plaquettes) && !any(a.syndromes.frontier)
        end
    end

    @testset "All small single X faults have syndrome-valid feedback routes" begin
        encoder = toric3d_slice_encoder(Toric3DSliceLattice(2,2,2))
        decoder = build_toric3d_slice_decoder(encoder.code)
        Hq = toric3d_slice_check_matrix(encoder.code)
        clean = sample_toric3d_faults(MersenneTwister(1),encoder)
        for (gate,op) in enumerate(encoder.operations), q in op.qubits
            record = with_toric3d_fault(encoder,clean,gate,q,:X)
            result = replay_toric3d_slice(encoder,record,decoder;record_history=true)
            @test !any(result.frame.z)
            for h in result.history
                @test slice_parity(Hq,h.recovery) == h.syndrome
                @test h.residual_syndrome == h.syndrome .⊻ slice_parity(Hq,h.route)
            end
        end
    end
end
