const _VERTEX_Z_AUDIT_FORMAT="h5-z-boundary-ambiguity-v1"
const _VERTEX_Z_AUDIT_FILENAMES=(
    json="z-boundary-audit.json",
    groups_csv="z-boundary-groups.csv",
    report="z-boundary-report.md",
)

_vertex_z_protocol_index(protocol)=something(findfirst(==(protocol),TORIC3D_VERTEX_PROTOCOLS),typemax(Int))
_vertex_z_class_text(cls)=string(Int(cls[1]),Int(cls[2]))

function _vertex_z_sorted_audits(audits)
    result=Toric3DVertexZAmbiguityAudit[audits...]
    isempty(result) && throw(ArgumentError("at least one Z-boundary audit is required"))
    identities=Set{Tuple{NTuple{3,Int},Symbol,Bool}}()
    for audit in result
        audit.protocol in TORIC3D_VERTEX_PROTOCOLS ||
            throw(ArgumentError("audit protocol is not canonical"))
        identity=(audit.dimensions,audit.protocol,audit.final_layer_only)
        identity in identities && throw(ArgumentError("duplicate Z-boundary audit identity"))
        push!(identities,identity)
        audit.mechanism_count>0 || throw(ArgumentError("Z-boundary audit is empty"))
        audit.mechanism_count==length(audit.signatures) ||
            throw(ArgumentError("Z-boundary audit mechanism count is inconsistent"))
        audit.group_count==length(audit.groups) ||
            throw(ArgumentError("Z-boundary audit group count is inconsistent"))
    end
    sort!(result;by=a->(a.dimensions...,_vertex_z_protocol_index(a.protocol),a.final_layer_only))
end

function _vertex_z_parity_payload(summary::Toric3DVertexZParitySummary)
    Dict{String,Any}(
        "mechanism_count"=>summary.mechanism_count,
        "group_count"=>summary.group_count,
        "mixed_group_count"=>summary.mixed_group_count,
        "unavoidable_count"=>summary.unavoidable_count,
        "baseline_failure_count"=>summary.baseline_failure_count,
        "control_failure_count"=>summary.unavoidable_count,
    )
end

function _vertex_z_group_payload(group::Toric3DVertexZObservableGroup)
    Dict{String,Any}(
        "syndrome_indices"=>collect(group.syndrome_indices),
        "upper_parity_odd"=>group.upper_parity_odd,
        "class_counts"=>collect(group.class_counts),
        "selected_class"=>collect(group.selected_class),
        "baseline_class"=>collect(group.baseline_class),
        "mixed"=>_vertex_z_mixed(group.class_counts),
        "unavoidable_count"=>group.unavoidable_count,
        "baseline_failure_count"=>group.baseline_failure_count,
        "control_failure_count"=>group.unavoidable_count,
    )
end

function _vertex_z_audit_dict(audit::Toric3DVertexZAmbiguityAudit)
    groups=sort!(collect(audit.groups);by=g->(g.upper_parity_odd,g.syndrome_indices))
    group_payloads=_vertex_z_group_payload.(groups)
    t,u,v=audit.dimensions
    Dict{String,Any}(
        "dimensions"=>[t,u,v],
        "protocol"=>string(audit.protocol),
        "final_layer_only"=>audit.final_layer_only,
        "encoder_operation_count"=>6*t*u*v,
        "active_output_fault_count"=>audit.mechanism_count,
        "mechanism_count"=>audit.mechanism_count,
        "group_count"=>audit.group_count,
        "mixed_group_count"=>audit.mixed_group_count,
        "unavoidable_count"=>audit.unavoidable_count,
        "unavoidable_rate"=>audit.unavoidable_count/audit.mechanism_count,
        "baseline_failure_count"=>audit.baseline_failure_count,
        "control_failure_count"=>audit.unavoidable_count,
        "parity"=>Dict(
            "even"=>_vertex_z_parity_payload(audit.parity.even),
            "odd"=>_vertex_z_parity_payload(audit.parity.odd),
        ),
        "groups"=>group_payloads,
        "mixed_groups"=>filter(g->g["mixed"],group_payloads),
    )
end

"""Create the deterministic machine-readable H+5 Z-boundary ambiguity certificate."""
function toric3d_vertex_z_audit_payload(audits;provenance=Dict())
    ordered=_vertex_z_sorted_audits(audits)
    payload=Dict{String,Any}(
        "format"=>_VERTEX_Z_AUDIT_FORMAT,
        "status"=>"diagnostic_only",
        "single_fault_prior"=>"equal single-location Z-fault priors",
        "decoder_observables"=>["final star syndrome","observed upper parity"],
        "logical_z_label"=>"state-independent residual logical Z; not |00_L> infidelity",
        "threshold_eligible"=>false,
        "protocols"=>string.(collect(TORIC3D_VERTEX_PROTOCOLS)),
        "source_hashes"=>_vertex_source_hashes(),
        "backend_versions"=>_vertex_versions(),
        "audits"=>_vertex_z_audit_dict.(ordered),
    )
    extra=Dict{String,Any}(string(k)=>v for (k,v) in pairs(provenance))
    isempty(intersect(Set(keys(payload)),Set(keys(extra)))) ||
        throw(ArgumentError("provenance cannot overwrite scientific keys"))
    merge!(payload,extra)
end

function _vertex_z_write_groups_csv(path,audits)
    _toric3d_atomic_write(path) do io
        println(io,"T,Lu,Lv,protocol,final_layer_only,upper_parity_odd,syndrome_indices,class_00,class_01,class_10,class_11,selected_class,baseline_class,mixed,unavoidable_count,baseline_failure_count,control_failure_count")
        for audit in audits
            groups=sort!(collect(audit.groups);by=g->(g.upper_parity_odd,g.syndrome_indices))
            for group in groups
                values=(audit.dimensions...,audit.protocol,audit.final_layer_only,
                    group.upper_parity_odd,join(group.syndrome_indices,';'),
                    group.class_counts...,_vertex_z_class_text(group.selected_class),
                    _vertex_z_class_text(group.baseline_class),_vertex_z_mixed(group.class_counts),
                    group.unavoidable_count,group.baseline_failure_count,group.unavoidable_count)
                println(io,join(_toric3d_csv.(values),','))
            end
        end
    end
end

function _vertex_z_write_report(path,audits)
    _toric3d_atomic_write(path) do io
        println(io,"# H+5 Z-boundary ambiguity diagnostic")
        println(io)
        println(io,"Status: `diagnostic_only`.")
        println(io)
        println(io,"Logical Z is a state-independent diagnostic, not `|00_L>` infidelity. The decoder-visible information is the final star syndrome and observed upper parity.")
        println(io)
        println(io,"| T | Lu | Lv | protocol | scope | locations | mixed groups | avoidable baseline failures | control / unavoidable floor |")
        println(io,"|---:|---:|---:|:---|:---|---:|---:|---:|---:|")
        for audit in audits
            scope=audit.final_layer_only ? "final layer" : "all active outputs"
            avoidable=audit.baseline_failure_count-audit.unavoidable_count
            println(io,"| $(audit.dimensions[1]) | $(audit.dimensions[2]) | $(audit.dimensions[3]) | `$(audit.protocol)` | $scope | $(audit.mechanism_count) | $(audit.mixed_group_count) | $avoidable | $(audit.unavoidable_count) |")
        end
        println(io)
        println(io,"The circuit-aware control removes the avoidable baseline mismatch but retains the information-theoretically unavoidable floor. It cannot establish a deterministic threshold for the current boundary instrument.")
        println(io)
        println(io,"Counts assume equal single-location Z-fault priors; tied logical classes use deterministic order `00, 01, 10, 11`.")
    end
end

"""Atomically write JSON, group CSV, and Markdown views of the ambiguity certificate."""
function write_toric3d_vertex_z_ambiguity_audit(audits,dir::AbstractString;provenance=Dict())
    ordered=_vertex_z_sorted_audits(audits)
    payload=toric3d_vertex_z_audit_payload(ordered;provenance)
    root=abspath(dir);mkpath(root)
    paths=(json=joinpath(root,_VERTEX_Z_AUDIT_FILENAMES.json),
        groups_csv=joinpath(root,_VERTEX_Z_AUDIT_FILENAMES.groups_csv),
        report=joinpath(root,_VERTEX_Z_AUDIT_FILENAMES.report))
    _vertex_write_json(paths.json,payload)
    _vertex_z_write_groups_csv(paths.groups_csv,ordered)
    _vertex_z_write_report(paths.report,ordered)
    paths
end
