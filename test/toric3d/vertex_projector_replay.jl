using Test,Random,TopoNoise
function h5_test_outcome(out,c)
    @test !any(TopoNoise._vertex_syndrome(star_checks(c),out.residual.z))
    @test !any(TopoNoise._vertex_syndrome(plaquette_checks(c),out.residual.x))
    @test out.logical_x==Tuple(TopoNoise._vertex_syndrome(toric3d_encoded_logical_z(c),out.residual.x))
    @test out.logical_z==Tuple(TopoNoise._vertex_syndrome(toric3d_encoded_logical_x(c),out.residual.z))
    @test out.pauli_class==sum(Int(b)<<(i-1) for (i,b) in enumerate((out.logical_x...,out.logical_z...)))
    @test out.upper_parity_odd==isodd(count(out.boundary.observed))
    for h in out.history
        @test TopoNoise._vertex_syndrome(h.window.matrix,h.full_recovery)==h.syndrome
        @test h.post_syndrome==h.syndrome.⊻TopoNoise._vertex_syndrome(h.window.matrix,h.route)
        @test 1<=h.checkpoint<c.T
        @test all(q->data_qubit_coordinate(c,q)[2]<=h.checkpoint,h.window.qubits)
    end
end
@testset "H+5 pre-final replay exposes the exact boundary-completed state" begin
    c=Toric3DVertexProjectorCode(3,2,3);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.1)
    clean=sample_toric3d_vertex_projector_record(MersenneTwister(41),e)
    step=e.checkpoints[end-1]+1
    fault=with_toric3d_vertex_projector_fault(clean,e,step,e.operations[step].qubits[1],:Z)
    for r in (clean,fault)
        pre=replay_toric3d_vertex_projector_pre_final(e,r,f;
            protocol=:window_complete,channel=:z_only)
        out=replay_toric3d_vertex_projector(e,r,f;
            protocol=:window_complete,channel=:z_only)
        @test pre.upper_parity_odd==isodd(count(pre.boundary.observed))
        @test pre.frame.z==out.boundary.z
        @test pre.before_boundary==out.before_boundary
        @test [(h.checkpoint,h.channel,h.signature) for h in pre.history]==
              [(h.checkpoint,h.channel,h.signature) for h in out.history]
        @test out.upper_parity_odd==pre.upper_parity_odd
    end
end
@testset "H+5 four-arm replay shares one immutable fault and boundary record" begin
    @test isdefined(TopoNoise,:replay_toric3d_vertex_projector)
    if isdefined(TopoNoise,:replay_toric3d_vertex_projector)
        @test TORIC3D_VERTEX_PROTOCOLS==(:none,:one_step,:window_complete,:full_prefix)
        @test TORIC3D_VERTEX_REPLAY_CHANNELS==(:x_only,:z_only,:both)
        c=Toric3DVertexProjectorCode(3,2,3);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.1)
        for channel in (:x_only,:z_only,:both)
            r=sample_toric3d_vertex_projector_record(MersenneTwister(94),e;p_x=channel==:z_only ? 0. : .1,p_z=channel==:x_only ? 0. : .1)
            saved=(copy(r.x),copy(r.z),copy(r.upper_branch))
            for protocol in TORIC3D_VERTEX_PROTOCOLS
                out=replay_toric3d_vertex_projector(e,r,f;protocol,channel)
                @test out.boundary.clean_branch==r.upper_branch
                @test length(out.history)==(protocol==:none ? 0 : (channel==:both ? 4 : 2))
                @test (r.x,r.z,r.upper_branch)==saved
                @test (channel==:x_only ? out.state_failure==any(out.logical_x) : isnothing(out.state_failure))
                h5_test_outcome(out,c)
            end
        end
        r=sample_toric3d_vertex_projector_record(MersenneTwister(1),e;p_x=1,p_z=1)
        @test_throws ArgumentError replay_toric3d_vertex_projector(e,r,f;protocol=:none,channel=:invalid)
        @test_throws ArgumentError replay_toric3d_vertex_projector(e,r,f;protocol=:none,channel=:x_only)
        @test_throws ArgumentError replay_toric3d_vertex_projector(e,r,f;protocol=:none,channel=:z_only)
        @test_throws ArgumentError replay_toric3d_vertex_projector(e,r,f;protocol=:invalid,channel=:both)
    end
end
if isdefined(TopoNoise,:replay_toric3d_vertex_projector)
@testset "H+5 exhaustive single-output faults have validated histories" begin
    c=Toric3DVertexProjectorCode(2,2,2);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.01)
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(2),e)
    for (step,op) in enumerate(e.operations),q in op.qubits,pauli in (:X,:Z)
        r=with_toric3d_vertex_projector_fault(zero,e,step,q,pauli)
        for protocol in TORIC3D_VERTEX_PROTOCOLS
            out=replay_toric3d_vertex_projector(e,r,f;protocol,channel=pauli==:X ? :x_only : :z_only)
            h5_test_outcome(out,c)
        end
    end
end
@testset "H+5 zero noise prepares the target state in both channels" begin
    for dims in ((1,2,2),(2,2,2),(3,2,3))
        c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=0)
        @test verify_toric3d_vertex_projector_encoder(e)
        for seed in 1:8,channel in TORIC3D_VERTEX_REPLAY_CHANNELS,protocol in TORIC3D_VERTEX_PROTOCOLS
            r=sample_toric3d_vertex_projector_record(MersenneTwister(seed),e)
            out=replay_toric3d_vertex_projector(e,r,f;protocol,channel)
            @test out.pauli_class==0
            @test !out.logical_x_failure && !out.logical_z_failure
            @test !any(out.residual.x)&&!any(out.residual.z)
            c.T==1 && @test isempty(out.history)
        end
    end
end
@testset "H+5 final recovery reports state-independent logical Pauli classes" begin
    c=Toric3DVertexProjectorCode(2,2,2);f=build_toric3d_vertex_decoders(c;p=.1);n=data_qubit_count(c)
    for mask in 0:15
        frame=PauliFrame(falses(n),falses(n))
        for i in 1:2
            isodd(mask>>(i-1)) && (frame.x[toric3d_encoded_logical_x(c)[i]].⊻=true)
            isodd(mask>>(i+1)) && (frame.z[toric3d_encoded_logical_z(c)[i]].⊻=true)
        end
        result=TopoNoise._vertex_final_recovery(frame,f)
        @test result.pauli_class==mask
        @test !any(result.x_recovery)&&!any(result.z_decision.full_recovery)
    end
end
@testset "H+5 prepared-state verification is distinct from syndrome-valid decoding" begin
    c=Toric3DVertexProjectorCode(2,2,2);f=build_toric3d_vertex_decoders(c;p=.1);n=data_qubit_count(c)
    x=falses(n);x[toric3d_encoded_logical_x(c)[1]].=true
    result=TopoNoise._vertex_final_recovery(PauliFrame(x,falses(n)),f)
    @test !any(result.x_syndrome)
    @test !any(result.x_recovery)
    @test result.logical_x==(true,false)
    # For :both every recorded pair of decisions sees the same pre-correction frame.
    e=toric3d_vertex_projector_encoder(c);r=sample_toric3d_vertex_projector_record(MersenneTwister(9),e;p_x=.15,p_z=.15)
    for protocol in (:one_step,:window_complete)
        out=replay_toric3d_vertex_projector(e,r,f;protocol,channel=:both)
        frame=PauliFrame(falses(n),falses(n));first=1
        for t in 1:c.T-1
            propagate_toric3d_vertex_projector_segment!(frame,e,r,first,e.checkpoints[t]);first=e.checkpoints[t]+1
            pair=filter(h->h.checkpoint==t,out.history)
            @test length(pair)==2
            for h in pair
                @test h.syndrome==TopoNoise._vertex_syndrome(h.window.checks,h.channel==:x_only ? frame.x : frame.z)
            end
            for h in reverse(pair)
                bits=h.channel==:x_only ? frame.x : frame.z;bits[h.window.qubits].⊻=h.route
            end
        end
        propagate_toric3d_vertex_projector_segment!(frame,e,r,first,length(e.operations))
        @test frame==out.before_boundary
    end
end
end
