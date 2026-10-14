using Test
using TopoNoise

function h5_z_control_audit_fixture(audit;dimensions=audit.dimensions,protocol=audit.protocol)
    source=first(filter(g->g.upper_parity_odd,audit.groups))
    group=Toric3DVertexZObservableGroup(source.syndrome_indices,true,
        (1,3,3,0),(false,true),source.baseline_class,4,
        source.baseline_failure_count)
    even=Toric3DVertexZParitySummary(0,0,0,0,0)
    odd=Toric3DVertexZParitySummary(7,1,1,4,group.baseline_failure_count)
    Toric3DVertexZAmbiguityAudit(dimensions,protocol,true,
        Toric3DVertexZFaultSignature[],[group],7,1,1,4,
        group.baseline_failure_count,(even=even,odd=odd))
end

@testset "H+5 circuit-aware Z control reaches the audited ambiguity floor" begin
    expected=Dict(3=>(8,34),5=>(14,131),7=>(20,276))
    for L in (3,5,7)
        c=Toric3DVertexProjectorCode(L,L,L)
        e=toric3d_vertex_projector_encoder(c)
        f=build_toric3d_vertex_decoders(c;p=.01)
        audit=audit_toric3d_vertex_z_boundary(e,f;
            protocol=:none,final_layer_only=true)
        decoder=build_toric3d_vertex_circuit_z_decoder(f,audit)
        zero=Toric3DVertexProjectorFaultRecord(c.T,c.Lu,c.Lv,
            falses(2,length(e.operations)),falses(2,length(e.operations)),
            falses(c.Lu*c.Lv))
        outcomes=Toric3DVertexZControlOutcome[]
        for signature in audit.signatures
            record=with_toric3d_vertex_projector_fault(zero,e,signature.step,
                signature.qubit,:Z)
            pre=replay_toric3d_vertex_projector_pre_final(e,record,f;
                protocol=:none,channel=:z_only,record_history=false)
            push!(outcomes,finalize_toric3d_vertex_z_control(pre,decoder))
        end
        floor,baseline=expected[L]
        @test count(out->out.control_logical_z_failure,outcomes)==audit.unavoidable_count
        @test audit.unavoidable_count==floor
        @test count(out->out.baseline_logical_z_failure,outcomes)==baseline
        @test all(out->!any(out.residual_syndrome),outcomes)
        @test all(out->!out.decision.fallback,outcomes)
    end
end

@testset "H+5 parity-erasure estimates retain acceptance cost and uncertainty" begin
    parity=Bool[false,false,true,true,true]
    failures=Bool[false,true,true,false,true]
    estimate=toric3d_vertex_parity_erasure_estimate(parity,failures)
    @test estimate.shots==5
    @test estimate.accepted==2
    @test estimate.rejected==3
    @test estimate.failures==3
    @test estimate.accepted_failures==1
    @test estimate.acceptance_rate==2/5
    @test estimate.acceptance_ci==TopoNoise._toric3d_wilson(2,5)
    @test estimate.conditional_rate==1/2
    @test estimate.conditional_ci==TopoNoise._toric3d_wilson(1,2)

    all_odd=toric3d_vertex_parity_erasure_estimate(trues(4),Bool[true,false,true,false])
    @test all_odd.accepted==0
    @test all_odd.rejected==4
    @test all_odd.conditional_rate===nothing
    @test all_odd.conditional_ci===nothing
    @test all(isfinite,all_odd.acceptance_ci)
    @test_throws DimensionMismatch toric3d_vertex_parity_erasure_estimate(
        Bool[false],Bool[false,true])
    @test_throws ArgumentError toric3d_vertex_parity_erasure_estimate(Bool[],Bool[])
end

@testset "H+5 circuit-aware Z control selects a deterministic logical class" begin
    c=Toric3DVertexProjectorCode(2,2,2)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)
    real=audit_toric3d_vertex_z_boundary(e,f;
        protocol=:none,final_layer_only=true)
    audit=h5_z_control_audit_fixture(real)
    decoder=build_toric3d_vertex_circuit_z_decoder(f,audit)
    group=only(audit.groups)
    syndrome=falses(length(star_checks(c)))
    syndrome[collect(group.syndrome_indices)].=true

    decision=decode_toric3d_vertex_circuit_z(decoder,syndrome,true)
    @test decision.selected_class==(false,true)
    @test decision.class_counts==(1,3,3,0)
    @test !decision.fallback
    @test TopoNoise._vertex_syndrome(star_checks(c),decision.full_recovery)==syndrome

    unseen=falses(length(syndrome))
    fallback=decode_toric3d_vertex_circuit_z(decoder,unseen,false)
    @test fallback.fallback
    @test fallback.full_recovery==fallback.base.full_recovery

    @test Set(keys(decoder.logical_adjustments))==
        Set(((false,false),(false,true),(true,false),(true,true)))
    for adjustment in values(decoder.logical_adjustments)
        @test !any(TopoNoise._vertex_syndrome(star_checks(c),adjustment))
    end

    @test_throws DimensionMismatch build_toric3d_vertex_circuit_z_decoder(f,
        h5_z_control_audit_fixture(real;dimensions=(3,2,2)))
    @test_throws ArgumentError build_toric3d_vertex_circuit_z_decoder(f,
        h5_z_control_audit_fixture(real;protocol=:invalid))
end
