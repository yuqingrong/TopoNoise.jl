using TopoNoise

const RESULT_DIR=@__DIR__
const EXPECTED=Dict(
    3=>(locations=99,mixed=5,baseline=34,floor=8,even=27,odd=72),
    5=>(locations=275,mixed=9,baseline=131,floor=14,even=75,odd=200),
    7=>(locations=539,mixed=13,baseline=276,floor=20,even=147,odd=392),
)

signature_row(s)=(s.step,s.slot,s.qubit,s.final_layer,s.syndrome_indices,
    s.upper_parity_odd,s.logical_z)

function control_failures(e,f,audit)
    decoder=build_toric3d_vertex_circuit_z_decoder(f,audit)
    c=e.code
    zero=Toric3DVertexProjectorFaultRecord(c.T,c.Lu,c.Lv,
        falses(2,length(e.operations)),falses(2,length(e.operations)),
        falses(c.Lu*c.Lv))
    failures=0
    for signature in audit.signatures
        record=with_toric3d_vertex_projector_fault(zero,e,signature.step,
            signature.qubit,:Z)
        pre=replay_toric3d_vertex_projector_pre_final(e,record,f;
            protocol=audit.protocol,channel=:z_only,record_history=false)
        outcome=finalize_toric3d_vertex_z_control(pre,decoder)
        !any(outcome.residual_syndrome) || error("control left final syndrome")
        failures+=outcome.control_logical_z_failure
    end
    failures
end

for key in ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS")
    get(ENV,key,"")=="1" || error("$key must equal 1")
end
get(ENV,"JULIA_CONDAPKG_BACKEND","")=="Null" ||
    error("JULIA_CONDAPKG_BACKEND must equal Null")

audits=Toric3DVertexZAmbiguityAudit[]
for L in (3,5,7)
    c=Toric3DVertexProjectorCode(L,L,L)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)
    by_protocol=Dict(protocol=>audit_toric3d_vertex_z_boundary(e,f;
        protocol,final_layer_only=true) for protocol in TORIC3D_VERTEX_PROTOCOLS)
    reference=signature_row.(by_protocol[:none].signatures)
    expected=EXPECTED[L]
    for protocol in TORIC3D_VERTEX_PROTOCOLS
        audit=by_protocol[protocol]
        signature_row.(audit.signatures)==reference ||
            error("final-layer signatures differ for L=$L protocol=$protocol")
        audit.mechanism_count==expected.locations || error("location count changed")
        audit.mixed_group_count==expected.mixed || error("mixed-group count changed")
        audit.baseline_failure_count==expected.baseline || error("baseline count changed")
        audit.unavoidable_count==expected.floor || error("unavoidable floor changed")
        audit.parity.even.mechanism_count==expected.even || error("even count changed")
        audit.parity.odd.mechanism_count==expected.odd || error("odd count changed")
        audit.parity.even.mixed_group_count==0 || error("even sector became ambiguous")
        control_failures(e,f,audit)==expected.floor || error("control missed audited floor")
        push!(audits,audit)
        println("L=$L protocol=$protocol locations=$(audit.mechanism_count) mixed=$(audit.mixed_group_count) baseline=$(audit.baseline_failure_count) control=$(audit.unavoidable_count)")
    end
end

implementation_head=readchomp(`git rev-parse HEAD`)
provenance=Dict{String,Any}(
    "spec_commit"=>"cb2ec97",
    "implementation_head"=>implementation_head,
    "execution_date"=>"2026-10-14",
    "single_thread_environment"=>Dict(key=>ENV[key] for key in
        ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS",
         "VECLIB_MAXIMUM_THREADS","JULIA_CONDAPKG_BACKEND","JULIA_PYTHONCALL_EXE")),
    "single_fault_prior_detail"=>"each active output has the same first-order Z-fault weight",
    "preliminary_result"=>"../toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7",
    "historical_full_scan"=>"../toric3d-h5-diagonal-feedback-2026-10-13/production (preserved; not resumed)",
    "note_tex_sha256"=>"3a40710c2ab6bd43e53f5f7703b7c8716b209dc1e060a5092f17e0edce073f3a",
)
write_toric3d_vertex_z_ambiguity_audit(audits,RESULT_DIR;provenance)

runtime=Dict{String,Any}(
    "julia_version"=>string(VERSION),
    "python_version"=>TopoNoise.PythonCall.pyconvert(String,
        TopoNoise.PythonCall.pyimport("sys").version),
    "backend_versions"=>TopoNoise._vertex_versions(),
    "implementation_head"=>implementation_head,
    "spec_commit"=>"cb2ec97",
    "project"=>something(Base.active_project(),""),
    "kernel"=>string(Sys.KERNEL),
    "architecture"=>string(Sys.ARCH),
    "environment"=>provenance["single_thread_environment"],
)
TopoNoise._vertex_write_json(joinpath(RESULT_DIR,"runtime-environment.json"),runtime)
