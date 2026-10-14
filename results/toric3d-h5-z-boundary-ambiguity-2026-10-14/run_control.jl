using Random
using TopoNoise

const RESULT_DIR=@__DIR__
const MASTER_SEED=1234
const SHOTS=1000
const RATES=(1e-4,1e-3)

csv_value(::Nothing)=""
csv_value(value)=TopoNoise._toric3d_csv(value)
interval_values(::Nothing)=(nothing,nothing)
interval_values(value)=value

for key in ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS")
    get(ENV,key,"")=="1" || error("$key must equal 1")
end
get(ENV,"JULIA_CONDAPKG_BACKEND","")=="Null" ||
    error("JULIA_CONDAPKG_BACKEND must equal Null")

master=Xoshiro(MASTER_SEED)
point_seeds=Dict((L,p)=>rand(master,UInt64) for L in (3,5,7) for p in RATES)
rows=Vector{Vector{Any}}()
for L in (3,5,7)
    c=Toric3DVertexProjectorCode(L,L,L)
    e=toric3d_vertex_projector_encoder(c)
    f=build_toric3d_vertex_decoders(c;p=.01)
    decoders=Dict(protocol=>build_toric3d_vertex_circuit_z_decoder(f,
        audit_toric3d_vertex_z_boundary(e,f;protocol,final_layer_only=true))
        for protocol in TORIC3D_VERTEX_PROTOCOLS)
    for p in RATES
        point_seed=point_seeds[(L,p)]
        rng=Xoshiro(point_seed)
        parity=Dict(protocol=>Bool[] for protocol in TORIC3D_VERTEX_PROTOCOLS)
        baseline=Dict(protocol=>Bool[] for protocol in TORIC3D_VERTEX_PROTOCOLS)
        control=Dict(protocol=>Bool[] for protocol in TORIC3D_VERTEX_PROTOCOLS)
        fallbacks=Dict(protocol=>0 for protocol in TORIC3D_VERTEX_PROTOCOLS)
        for _ in 1:SHOTS
            record=sample_toric3d_vertex_projector_record(rng,e;p_x=0,p_z=p)
            for protocol in TORIC3D_VERTEX_PROTOCOLS
                pre=replay_toric3d_vertex_projector_pre_final(e,record,f;
                    protocol,channel=:z_only,record_history=false)
                outcome=finalize_toric3d_vertex_z_control(pre,decoders[protocol])
                push!(parity[protocol],pre.upper_parity_odd)
                push!(baseline[protocol],outcome.baseline_logical_z_failure)
                push!(control[protocol],outcome.control_logical_z_failure)
                fallbacks[protocol]+=outcome.decision.fallback
            end
        end
        for protocol in TORIC3D_VERTEX_PROTOCOLS
            base_est=toric3d_vertex_parity_erasure_estimate(parity[protocol],baseline[protocol])
            control_est=toric3d_vertex_parity_erasure_estimate(parity[protocol],control[protocol])
            base_ci=TopoNoise._toric3d_wilson(base_est.failures,SHOTS)
            control_ci=TopoNoise._toric3d_wilson(control_est.failures,SHOTS)
            base_conditional=interval_values(base_est.conditional_ci)
            control_conditional=interval_values(control_est.conditional_ci)
            push!(rows,Any[
                "parity-erasure diagnostic",L,L,L,p,SHOTS,MASTER_SEED,string(point_seed),
                protocol,base_est.accepted,base_est.rejected,base_est.failures,
                control_est.failures,base_est.accepted_failures,
                control_est.accepted_failures,base_est.acceptance_rate,
                base_est.acceptance_ci...,base_est.failures/SHOTS,base_ci...,
                control_est.failures/SHOTS,control_ci...,
                base_est.conditional_rate,base_conditional...,
                control_est.conditional_rate,control_conditional...,
                fallbacks[protocol],
            ])
        end
        println("L=$L p=$p seed=$point_seed completed $SHOTS paired records")
    end
end

path=joinpath(RESULT_DIR,"parity-erasure-control.csv")
header=["label","T","Lu","Lv","p","shots","master_seed","point_seed","protocol",
    "accepted_even","rejected_odd","baseline_failures","control_failures",
    "accepted_baseline_failures","accepted_control_failures","acceptance_rate",
    "acceptance_ci_low","acceptance_ci_high","baseline_rate","baseline_ci_low",
    "baseline_ci_high","control_rate","control_ci_low","control_ci_high",
    "conditional_baseline_rate","conditional_baseline_ci_low","conditional_baseline_ci_high",
    "conditional_control_rate","conditional_control_ci_low","conditional_control_ci_high",
    "single_fault_table_fallbacks"]
all(length(row)==length(header) for row in rows) || error("control CSV schema mismatch")
TopoNoise._toric3d_atomic_write(path) do io
    println(io,join(header,','))
    for row in rows
        println(io,join(csv_value.(row),','))
    end
end
