using Test, Random, TopoNoise
if !isdefined(@__MODULE__, :h5_physical_tables)
    include("vertex_projector_test_helpers.jl")
end
# This reference uses matrix-conjugation Pauli tables, direct support parities,
# and applies the returned physical corrections; it never calls frame propagation.
h5_parities(checks,labels,bit)=BitVector([isodd(sum((labels[q]&bit)!=0 for q in support)) for support in checks])
function h5_reference_history(e,r,f,protocol,channel,tables)
    out=replay_toric3d_vertex_projector(e,r,f;protocol,channel)
    labels=zeros(Int,data_qubit_count(e.code));first=1;m=e.code.Lu*e.code.Lv
    for t in 1:e.code.T
        h5_physical_segment!(labels,e,r,first,e.checkpoints[t],tables);first=e.checkpoints[t]+1
        for h in filter(h->h.checkpoint==t,out.history)
            @testset "checkpoint=$t channel=$(h.channel) support=$(h.window.qubits) syndrome=$(h.syndrome) route=$(findall(h.route)) branch=$(collect(r.upper_branch))" begin
                bit=h.channel==:x_only ? 1 : 2
                @test h5_parities(h.window.checks,labels,bit)==h.syndrome
                @test all(q->data_qubit_coordinate(e.code,q)[2]<=t,h.window.qubits)
                corrected=zeros(Int,length(labels));corrected[h.window.qubits[findall(h.full_recovery)]].=bit
                @test h5_parities(h.window.checks,corrected,bit)==h.syndrome
                labels[h.window.qubits[findall(h.route)]].⊻=bit
                @test h5_parities(h.window.checks,labels,bit)==h.post_syndrome
                if protocol!=:one_step
                    prepared=h.channel==:x_only ? plaquette_checks(e.code)[1:3m*t] : star_checks(e.code)[1:(t+1)*m]
                    @test !any(h5_parities(prepared,labels,bit))
                elseif h.channel==:z_only
                    @test !any(h.post_syndrome[h.window.retired_rows])
                end
            end
        end
    end
    @test labels==Int.(out.before_boundary.x).+2Int.(out.before_boundary.z)
    a=BitVector(r.upper_branch).⊻h5_parities(f.boundary.upper_checks,labels,2)
    physical=falses(length(labels));clean=copy(physical)
    for i in eachindex(a)
        a[i] && (physical[f.boundary.paths[i]].⊻=true)
        r.upper_branch[i] && (clean[f.boundary.paths[i]].⊻=true)
    end
    @test out.boundary.observed==a
    @test out.boundary.physical_correction==physical
    @test out.boundary.clean_correction==clean
    labels[findall(physical)].⊻=2;labels[findall(clean)].⊻=2
    @test h5_parities(plaquette_checks(e.code),labels,1)==out.x_syndrome
    @test h5_parities(star_checks(e.code),labels,2)==out.z_syndrome
    labels[findall(out.x_recovery)].⊻=1;labels[findall(out.z_decision.full_recovery)].⊻=2
    @test !any(h5_parities(plaquette_checks(e.code),labels,1))
    @test !any(h5_parities(star_checks(e.code),labels,2))
    @test labels==Int.(out.residual.x).+2Int.(out.residual.z)
    @test Tuple(h5_parities(toric3d_encoded_logical_z(e.code),labels,1))==out.logical_x
    @test Tuple(h5_parities(toric3d_encoded_logical_x(e.code),labels,2))==out.logical_z
end
@testset "H+5 exhaustive single faults satisfy every checkpoint invariant" begin
    c=Toric3DVertexProjectorCode(2,2,2);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.01)
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(31),e);tables=h5_physical_tables()
    for (step,op) in enumerate(e.operations),q in op.qubits,pauli in (:X,:Z)
        @testset "gate=$step output=$q Pauli=$pauli" begin
            r=with_toric3d_vertex_projector_fault(zero,e,step,q,pauli)
            for protocol in TORIC3D_VERTEX_PROTOCOLS
                h5_reference_history(e,r,f,protocol,pauli==:X ? :x_only : :z_only,tables)
            end
        end
    end
end
@testset "H+5 randomized period-two histories match physical Pauli evolution" begin
    tables=h5_physical_tables()
    for dims in ((3,2,2),(3,2,3),(3,3,2)),seed in (17,91,483)
        c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c);f=build_toric3d_vertex_decoders(c;p=.15)
        for channel in TORIC3D_VERTEX_REPLAY_CHANNELS
            r=sample_toric3d_vertex_projector_record(MersenneTwister(seed),e;p_x=channel==:z_only ? 0. : .15,p_z=channel==:x_only ? 0. : .15)
            for protocol in TORIC3D_VERTEX_PROTOCOLS
                h5_reference_history(e,r,f,protocol,channel,tables)
            end
        end
    end
end
