using TopoNoise, Random, LinearAlgebra, TOML, Statistics
import PythonCall

const OUTPUT = "/private/tmp/h6-matched-x-pilot.toml"
const STARTED = time()
const LIMIT_SECONDS = 360.0

function replay_three(encoder, faults, slice, two)
    saved = deepcopy(faults)
    a = replay_toric3d_x_feedback(encoder, faults, slice, two.final;
        feedback=false, routing=:one_step)
    b = replay_toric3d_x_feedback(encoder, faults, slice, two.final;
        feedback=true, routing=:one_step)
    c = replay_toric3d_h6_two_slice_x(encoder, faults, two; feedback=true)
    @assert faults == saved
    for result in (a,b,c)
        @assert all(getproperty(result,k) isa Bool for k in (:logical_1,:logical_2,:failure))
        @assert result.failure == (result.logical_1 | result.logical_2)
        @assert !any(result.frame.z)
        @assert !any(TopoNoise._toric3d_matrix_syndrome(two.final.matrix, result.residual))
    end
    @assert !any(c.syndrome) && !any(c.recovery)
    return nothing
end

function save(data)
    data["elapsed_seconds"] = time()-STARTED
    open(OUTPUT*".tmp","w") do io
        TOML.print(io,data;sorted=true)
    end
    mv(OUTPUT*".tmp",OUTPUT;force=true)
end

function main()
    @assert Threads.nthreads() == BLAS.get_num_threads() == 1
    data = Dict{String,Any}("purpose"=>"Performance pilot only; all samples excluded from production",
        "canonical_project"=>dirname(Base.active_project()),"shots_per_point"=>3,
        "arm_order"=>["final_only","one_step","two_slice"],
        "julia_version"=>string(VERSION), "julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),
        "python_executable"=>PythonCall.pyconvert(String,PythonCall.pyimport("sys").executable),
        "points"=>Dict{String,Any}[],"complete"=>false)
    save(data)
    for L in (9,11,13,15)
        time()-STARTED > LIMIT_SECONDS && break
        geometry_seconds = @elapsed begin
            code=Toric3DEncodedSlice(L,L,L)
            encoder=toric3d_encoded_encoder(code)
            slice=build_toric3d_slice_decoder(code.lattice)
        end
        for p in (.003,.015,.05)
            time()-STARTED > LIMIT_SECONDS && break
            build_seconds = @elapsed two=build_toric3d_h6_two_slice_x_decoder(code;p)
            rng=Random.Xoshiro(UInt64(2026101700+100L+round(Int,1000p)))
            warmup_seconds = @elapsed begin
                faults=sample_toric3d_faults(rng,encoder;p_x=p,p_z=0.0)
                replay_three(encoder,faults,slice,two)
            end
            measured=Float64[]
            for _ in 1:3
                push!(measured,@elapsed begin
                    faults=sample_toric3d_faults(rng,encoder;p_x=p,p_z=0.0)
                    replay_three(encoder,faults,slice,two)
                end)
            end
            push!(data["points"],Dict{String,Any}("L"=>L,"p"=>p,
                "geometry_build_seconds"=>geometry_seconds,"decoder_build_seconds"=>build_seconds,
                "warmup_seconds"=>warmup_seconds,"shot_seconds"=>measured,
                "seconds_per_triple"=>mean(measured),"records"=>3,
                "projected_10000_shots_hours"=>10000mean(measured)/3600))
            save(data)
            println("L=",L," p=",p," build=",round(build_seconds;digits=3),
                " s; warmup=",round(warmup_seconds;digits=3),
                " s; mean triple=",round(mean(measured);digits=4),
                " s; 10000 triples=",round(10000mean(measured)/3600;digits=2)," h")
            flush(stdout)
        end
    end
    data["complete"] = length(data["points"]) == 12
    save(data)
end

main()
