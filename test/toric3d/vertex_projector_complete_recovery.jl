using Test,Random,TopoNoise
function h5_test_complete_history(e,r,family,policy,channel)
    frame=PauliFrame(falses(data_qubit_count(e.code)),falses(data_qubit_count(e.code)))
    ds=getproperty(family,Symbol(channel==:x_only ? "x_" : "z_",policy==:window_complete ? "window" : "prefix"))
    first=1;m=e.code.Lu*e.code.Lv
    for t in 1:e.code.T-1
        propagate_toric3d_vertex_projector_segment!(frame,e,r,first,e.checkpoints[t]);first=e.checkpoints[t]+1
        d=ds[t];w=d.window;bits=channel==:x_only ? frame.x : frame.z
        s=TopoNoise._vertex_syndrome(w.checks,bits)
        recovery=channel==:x_only ? decode_toric3d_vertex_bposd(d,s) : decode_toric3d_vertex_z_complete(d,s).full_recovery
        @test TopoNoise._vertex_syndrome(w.matrix,recovery)==s
        bits[w.qubits].⊻=recovery
        checks=channel==:x_only ? plaquette_checks(e.code)[1:3m*t] : star_checks(e.code)[1:(t+1)*m]
        @test !any(TopoNoise._vertex_syndrome(checks,bits))
        @test all(q->data_qubit_coordinate(e.code,q)[2]<=t,w.qubits)
    end
end
@testset "H+5 BP+OSD priors and every complete recovery are validated" begin
    @test isdefined(TopoNoise,:build_toric3d_vertex_decoders)
    if isdefined(TopoNoise,:build_toric3d_vertex_decoders)
        c=Toric3DVertexProjectorCode(2,2,2)
        for p in (0.,.02,1.)
            f=build_toric3d_vertex_decoders(c;p)
            for d in vcat(f.x_window,f.x_prefix,[f.final_x])
                @test d.p==p
                @test d.p_eff==clamp(p,1e-12,1-1e-12)
                @test d.version=="2.4.1"
                @test TopoNoise._toric3d_settings_metadata(d.settings)==(max_iter=50,bp_method="minimum_sum",ms_scaling_factor=1.,schedule="parallel",omp_thread_count=1,osd_method="OSD_CS",osd_order=2)
                s=falses(size(d.window.matrix,1));x=decode_toric3d_vertex_bposd(d,s)
                @test !any(x)
                x[1]=true;@test_throws TopoNoise.Toric3DVertexInvalidRecovery TopoNoise._vertex_validate_recovery(d.window,s,x)
            end
            for d in vcat(f.z_window,f.z_prefix,[f.final_z])
                s=falses(size(d.window.matrix,1));x=decode_toric3d_vertex_z_complete(d,s).full_recovery
                x[1]=true;@test_throws TopoNoise.Toric3DVertexInvalidRecovery TopoNoise._vertex_validate_recovery(d.window,s,x)
            end
        end
    end
end
if isdefined(TopoNoise,:build_toric3d_vertex_decoders)
@testset "H+5 window-complete X recovery clears every completed plaquette" begin
    c=Toric3DVertexProjectorCode(3,2,2);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.02)
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(1),e)
    for (step,op) in enumerate(e.operations),q in op.qubits
        r=with_toric3d_vertex_projector_fault(zero,e,step,q,:X)
        h5_test_complete_history(e,r,f,:window_complete,:x_only)
    end
    for dims in ((3,2,3),(3,3,2))
        c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.09);rng=MersenneTwister(72)
        for _ in 1:12
            r=sample_toric3d_vertex_projector_record(rng,e;p_x=.09)
            h5_test_complete_history(e,r,f,:window_complete,:x_only)
        end
    end
end
@testset "H+5 window-complete Z recovery clears every completed star sheet" begin
    c=Toric3DVertexProjectorCode(3,2,2);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.02)
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(1),e)
    for (step,op) in enumerate(e.operations),q in op.qubits
        r=with_toric3d_vertex_projector_fault(zero,e,step,q,:Z)
        h5_test_complete_history(e,r,f,:window_complete,:z_only)
    end
    for d in (f.z_window[1],f.z_prefix[2],f.final_z)
        s=falses(size(d.window.matrix,1));s[[1,2]].=true
        @test decode_toric3d_vertex_z_complete(d,s).signature==[(1,2,(6,10))]
    end
    f=build_toric3d_vertex_decoders(Toric3DVertexProjectorCode(3,3,3);p=.1)
    d=f.z_window[1];s=falses(18);s[[10,14]].=true
    a=decode_toric3d_vertex_z_complete(d,s)
    @test count(p->p.second==19,a.paths)==2
    @test count(a.full_recovery)==2
    s[14]=false;a=decode_toric3d_vertex_z_complete(d,s)
    @test count(a.full_recovery)==1
end
@testset "H+5 full-prefix recovery never accesses future data" begin
    for dims in ((3,2,2),(3,2,3),(3,3,2))
        c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.1);rng=MersenneTwister(8)
        for ch in (:x_only,:z_only),_ in 1:12
            r=sample_toric3d_vertex_projector_record(rng,e;p_x=ch==:x_only ? .1 : 0.,p_z=ch==:z_only ? .1 : 0.)
            h5_test_complete_history(e,r,f,:full_prefix,ch)
            h5_test_complete_history(e,r,f,:window_complete,ch)
        end
    end
end
@testset "H+5 combined-channel decisions share one pre-correction frame and commute" begin
    c=Toric3DVertexProjectorCode(2,2,2);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.1)
    r=sample_toric3d_vertex_projector_record(MersenneTwister(7),e;p_x=.1,p_z=.1)
    frame=PauliFrame(falses(36),falses(36));propagate_toric3d_vertex_projector_segment!(frame,e,r,1,24)
    xw=f.x_window[1].window;zw=f.z_window[1].window
    sx=TopoNoise._vertex_syndrome(xw.checks,frame.x);sz=TopoNoise._vertex_syndrome(zw.checks,frame.z)
    x=decode_toric3d_vertex_bposd(f.x_window[1],sx);z=decode_toric3d_vertex_z_complete(f.z_window[1],sz).full_recovery
    a=PauliFrame(frame.x,frame.z);b=PauliFrame(frame.x,frame.z)
    a.x[xw.qubits].⊻=x;a.z[zw.qubits].⊻=z
    b.z[zw.qubits].⊻=z;b.x[xw.qubits].⊻=x
    @test a==b
    @test !any(TopoNoise._vertex_syndrome(xw.checks,a.x))
    @test !any(TopoNoise._vertex_syndrome(zw.checks,a.z))
end
end
