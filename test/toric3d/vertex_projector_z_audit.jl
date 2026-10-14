using Test
using TopoNoise

@testset "H+5 Z-boundary audit enumerates and groups decoder-visible signatures" begin
    c=Toric3DVertexProjectorCode(3,3,3)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)

    audit=audit_toric3d_vertex_z_boundary(e,f;
        protocol=:none,final_layer_only=true)
    @test length(audit.signatures)==99
    @test length(audit.groups)==36
    @test audit.mechanism_count==99
    @test audit.group_count==36
    @test audit.mixed_group_count==5
    @test audit.unavoidable_count==8
    @test audit.baseline_failure_count==34
    @test audit.parity.even.mixed_group_count==0
    @test audit.parity.even.mechanism_count==27
    @test audit.parity.odd.mechanism_count==72

    all_signatures=toric3d_vertex_z_fault_signatures(e,f;
        protocol=:none,final_layer_only=false)
    expected=[(step,slot,q) for (step,op) in enumerate(e.operations)
              for (slot,q) in enumerate(op.qubits)]
    @test length(all_signatures)==11*c.T*c.Lu*c.Lv
    @test [(s.step,s.slot,s.qubit) for s in all_signatures]==expected
    @test all(g->g.selected_class in ((false,false),(false,true),(true,false),(true,true)),
        audit.groups)
    classes=((false,false),(false,true),(true,false),(true,true))
    @test all(g->g.selected_class==classes[findfirst(==(maximum(g.class_counts)),g.class_counts)],
        audit.groups)
    @test_throws ArgumentError audit_toric3d_vertex_z_boundary(e,f;
        protocol=:invalid,final_layer_only=true)
end

@testset "H+5 Z-boundary signatures are branch and protocol invariant" begin
    c=Toric3DVertexProjectorCode(3,3,3)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)
    m=c.Lu*c.Lv;nops=length(e.operations)
    branches=[falses(m),BitVector([i in (1,2) for i in 1:m]),
        BitVector([i in (m-1,m) for i in 1:m])]
    first_final=e.checkpoints[end-1]+1
    locations=((first_final,1),(e.checkpoints[end],length(e.operations[end].qubits)))
    for (step,slot) in locations
        observed=[]
        for branch in branches
            zero=Toric3DVertexProjectorFaultRecord(c.T,c.Lu,c.Lv,
                falses(2,nops),falses(2,nops),branch)
            record=with_toric3d_vertex_projector_fault(zero,e,step,
                e.operations[step].qubits[slot],:Z)
            pre=replay_toric3d_vertex_projector_pre_final(e,record,f;
                protocol=:none,channel=:z_only,record_history=false)
            syndrome=TopoNoise._vertex_syndrome(star_checks(c),pre.frame.z)
            logical=Tuple(TopoNoise._vertex_syndrome(toric3d_encoded_logical_x(c),pre.frame.z))
            push!(observed,(copy(pre.frame.z),Tuple(findall(syndrome)),pre.upper_parity_odd,logical))
        end
        @test all(==(first(observed)),observed)
    end

    signature_row(s)=(s.step,s.slot,s.qubit,s.final_layer,s.syndrome_indices,
        s.upper_parity_odd,s.logical_z)
    reference=signature_row.(toric3d_vertex_z_fault_signatures(e,f;
        protocol=:none,final_layer_only=true))
    for protocol in TORIC3D_VERTEX_PROTOCOLS
        @test signature_row.(toric3d_vertex_z_fault_signatures(e,f;
            protocol,final_layer_only=true))==reference
    end

    other=Toric3DVertexProjectorCode(2,2,2)
    other_family=build_toric3d_vertex_decoders(other;p=.01)
    @test_throws DimensionMismatch toric3d_vertex_z_fault_signatures(e,other_family;
        protocol=:none,final_layer_only=true)
end

@testset "H+5 final-layer ambiguity floors are pinned at L=3,5,7" begin
    expected=Dict(
        3=>(99,5,8,34,27,72),
        5=>(275,9,14,131,75,200),
        7=>(539,13,20,276,147,392),
    )
    for L in (3,5,7)
        c=Toric3DVertexProjectorCode(L,L,L)
        e=toric3d_vertex_projector_encoder(c)
        f=build_toric3d_vertex_decoders(c;p=.01)
        audit=audit_toric3d_vertex_z_boundary(e,f;
            protocol=:none,final_layer_only=true)
        locations,mixed,floor,baseline,even,odd=expected[L]
        @test audit.mechanism_count==locations
        @test audit.mixed_group_count==mixed
        @test audit.unavoidable_count==floor
        @test audit.baseline_failure_count==baseline
        @test audit.parity.even.mechanism_count==even
        @test audit.parity.odd.mechanism_count==odd
        @test audit.parity.even.mixed_group_count==0
    end
end
