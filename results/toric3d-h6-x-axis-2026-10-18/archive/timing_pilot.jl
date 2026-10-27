using Random, TOML
include("AxisH6.jl")
const A=AxisH6

function main()
    results=Dict{String,Any}[]
    for channel in ("x_only","z_only"), L in (3,7,11)
        rates=channel == "x_only" ? [0.003,0.015,0.055] : [0.006,0.02,0.035]
        previous=nothing
        for p in rates
            context=nothing
            build=@elapsed context=A.build_context(L,p,channel;previous)
            previous=context
            rng=Xoshiro(0x6a018 + L + (channel == "z_only" ? 100 : 0))
            warmup=@elapsed A.three_outcomes(context,A.sample_record(rng,context))
            failures=zeros(Int,3)
            elapsed=@elapsed for _ in 1:20
                outcomes=A.three_outcomes(context,A.sample_record(rng,context))
                failures .+= [result.failure for result in values(outcomes)]
            end
            row=Dict{String,Any}("channel"=>channel,"L"=>L,"p"=>p,"timed_records"=>20,
                "build_seconds"=>build,"warmup_seconds"=>warmup,"elapsed_seconds"=>elapsed,
                "seconds_per_record"=>elapsed/20,"failures"=>failures)
            push!(results,row)
            println(row); flush(stdout)
            path=joinpath(@__DIR__,"timing-pilot.toml")
            open(path*".tmp","w") do io
                TOML.print(io,Dict("excluded_from_production"=>true,"results"=>results);sorted=true)
            end
            mv(path*".tmp",path;force=true)
        end
    end
end
main()
