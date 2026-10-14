using Test
using TopoNoise

function h5_z_payload_is_finite(value)
    if value isa AbstractDict
        return all(h5_z_payload_is_finite,values(value))
    elseif value isa AbstractVector || value isa Tuple
        return all(h5_z_payload_is_finite,value)
    elseif value isa Real
        return isfinite(value)
    end
    true
end

@testset "H+5 Z ambiguity artifacts are deterministic and complete" begin
    c=Toric3DVertexProjectorCode(3,3,3)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)
    audits=[audit_toric3d_vertex_z_boundary(e,f;
        protocol,final_layer_only=true) for protocol in reverse(TORIC3D_VERTEX_PROTOCOLS)]
    payload=toric3d_vertex_z_audit_payload(audits;
        provenance=Dict("study"=>"test"))

    @test payload["format"]=="h5-z-boundary-ambiguity-v1"
    @test payload["status"]=="diagnostic_only"
    @test occursin("equal",lowercase(payload["single_fault_prior"]))
    @test payload["backend_versions"]==Dict("pymatching"=>"2.4.0","ldpc"=>"2.4.1")
    @test !isempty(payload["source_hashes"])
    @test [a["protocol"] for a in payload["audits"]]==collect(string.(TORIC3D_VERTEX_PROTOCOLS))
    @test all(a->a["mechanism_count"]==99,payload["audits"])
    @test all(a->a["parity"]["even"]["mechanism_count"]==27,payload["audits"])
    @test all(a->a["parity"]["odd"]["mechanism_count"]==72,payload["audits"])
    @test all(a->length(a["mixed_groups"])==a["mixed_group_count"],payload["audits"])
    @test all(a->a["control_failure_count"]==a["unavoidable_count"],payload["audits"])
    @test h5_z_payload_is_finite(payload)

    mktempdir() do dir
        paths=write_toric3d_vertex_z_ambiguity_audit(audits,dir;
            provenance=Dict("study"=>"test"))
        before=Dict(path=>read(path) for path in values(paths))
        again=write_toric3d_vertex_z_ambiguity_audit(audits,dir;
            provenance=Dict("study"=>"test"))
        @test paths==again
        @test all(path->read(path)==before[path],values(paths))

        rows=readlines(paths.groups_csv)
        @test startswith(first(rows),"T,Lu,Lv,protocol,final_layer_only,upper_parity_odd")
        protocol_index=Dict(string(p)=>i for (i,p) in enumerate(TORIC3D_VERTEX_PROTOCOLS))
        row_key(row)=begin
            fields=split(row,',')
            unquote(value)=startswith(value,'"') ? value[2:end-1] : value
            syndrome_text=unquote(fields[7])
            syndrome=isempty(syndrome_text) ? () : Tuple(parse.(Int,split(syndrome_text,';')))
            (parse(Int,fields[1]),parse(Int,fields[2]),parse(Int,fields[3]),
                protocol_index[unquote(fields[4])],parse(Bool,fields[6]),syndrome)
        end
        keys=row_key.(rows[2:end])
        @test keys==sort(keys)
        report=read(paths.report,String)
        @test occursin("diagnostic_only",report)
        @test occursin("avoidable baseline",lowercase(report))
        @test occursin("unavoidable",lowercase(report))
        @test occursin("not `|00_L>` infidelity",report)
        @test occursin("cannot establish a deterministic threshold",lowercase(report))
    end

    @test_throws ArgumentError toric3d_vertex_z_audit_payload(audits;
        provenance=Dict("format"=>"collision"))
    @test_throws ArgumentError toric3d_vertex_z_audit_payload(Toric3DVertexZAmbiguityAudit[])
    @test_throws ArgumentError toric3d_vertex_z_audit_payload([first(audits),first(audits)])
end
