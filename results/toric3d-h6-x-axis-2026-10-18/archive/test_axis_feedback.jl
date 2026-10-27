using Test, Random

# Test the archive in isolation, including while its implementation is absent.
module AxisFeedbackUnderTest
using TopoNoise, Random, PythonCall
const TN = TopoNoise
for name in ("axis_geometry.jl", "axis_feedback.jl")
    source = joinpath(@__DIR__, name)
    isfile(source) && include(source)
end
end
const AF = AxisFeedbackUnderTest

@testset "physical-x feedback public API" begin
    @test isdefined(AF, :AxisRecord)
    @test isdefined(AF, :build_context)
    @test isdefined(AF, :sample_record)
    @test isdefined(AF, :three_outcomes)
end

parity(checks, bits) = BitVector(isodd(sum(bits[q] for q in support)) for support in checks)
function cleanrecord(context)
    AF.AxisRecord(falses(2context.code.Ly*context.code.Lz),
        falses(2,length(context.encoder.operations)), falses(2,length(context.encoder.operations)))
end

# Deliberately independent of the production propagation helper.
function handsegment!(x,z,context,record,first,last)
    for index in first:last
        operation = context.encoder.operations[index]
        if operation.gate == :H
            q = only(operation.qubits)
            x[q],z[q] = z[q],x[q]
        elseif operation.gate == :CNOT
            c,t = operation.qubits
            x[t] = xor(x[t],x[c])
            z[c] = xor(z[c],z[t])
        else
            error("unrecognized operation")
        end
        for (slot,q) in enumerate(operation.qubits)
            x[q] = xor(x[q],record.x[slot,index])
            z[q] = xor(z[q],record.z[slot,index])
        end
    end
end

function audit_outcome(context,record,arm,outcome)
    code = context.code
    x,z = falses(AF.nqubits(code)),falses(AF.nqubits(code))
    z[AF.slice_qubits(code,0)] .= record.seed
    bits = context.channel == "x_only" ? x : z
    history = outcome.history
    if context.channel == "z_only"
        seed = first(history)
        @test seed.kind == :seed
        @test seed.checkpoint == 0
        @test seed.pre_syndrome == parity(seed.checks,z)
        @test seed.qubits == AF.slice_qubits(code,0)
        bits[seed.qubits] .⊻= seed.recovery
        @test !any(parity(seed.checks,z))
    end
    firstop = 1
    for r in 1:code.T
        entries = filter(entry -> entry.checkpoint == r,history)
        if context.channel == "x_only" && arm == :one_step
            entry = only(entries)
            @test entry.kind == :x_endpoint
            @test entry.pre_syndrome == parity(entry.checks,x)
            @test all(q -> q in AF.plane_edges(code,r-1),vcat(entry.checks...))
        end
        handsegment!(x,z,context,record,firstop,context.encoder.checkpoints[r])
        firstop = context.encoder.checkpoints[r]+1
        for entry in entries
            if !(context.channel == "x_only" && arm == :one_step)
                @test entry.pre_syndrome == parity(entry.checks,bits)
            end
            @test entry.support == entry.qubits[findall(entry.recovery)]
            bits[entry.qubits] .⊻= entry.recovery
            @test entry.post_syndrome == parity(entry.post_checks,bits)
            if arm == :two_slice
                prefix = AF.AxisCode(r,code.Ly,code.Lz)
                checks = context.channel == "x_only" ? AF.plaquette_checks(prefix) : AF.star_checks(prefix)
                allowed = vcat(AF.slice_qubits(code,r-1),AF.slice_qubits(code,r))
                expected = filter(support -> any(q -> q in allowed,support),checks)
                @test entry.checks == expected
                @test entry.qubits == allowed
                @test !any(parity(checks,bits))
                @test !any(entry.prefix_syndrome)
            elseif context.channel == "z_only" && arm == :one_step
                @test entry.kind == :z_push
                @test entry.recovery == entry.pre_syndrome
                @test entry.qubits == [AF.edge(code,:x,r,u,v) for u in 0:code.Ly-1 for v in 0:code.Lz-1]
                prefixstars = AF.star_checks(AF.AxisCode(r,code.Ly,code.Lz))
                @test entry.checks == prefixstars[(r-1)*code.Ly*code.Lz+1:r*code.Ly*code.Lz]
                @test !any(parity(prefixstars[1:r*code.Ly*code.Lz],z))
            elseif context.channel == "x_only" && arm == :one_step
                @test entry.qubits == AF.plane_edges(code,r)
                @test length(entry.recovery) == 2code.Ly*code.Lz
            end
        end
    end
    @test x == outcome.frame.x
    @test z == outcome.frame.z
    checks = context.channel == "x_only" ? AF.plaquette_checks(code) : AF.star_checks(code)
    logicals = context.channel == "x_only" ? AF.logical_z(code) : AF.logical_x(code)
    @test outcome.syndrome == parity(checks,bits)
    @test parity(checks,outcome.recovery) == outcome.syndrome
    @test outcome.residual == bits .⊻ outcome.recovery
    @test !any(parity(checks,outcome.residual))
    flags = parity(logicals,outcome.residual)
    @test (outcome.logical_1,outcome.logical_2) == Tuple(flags)
    @test outcome.failure == any(flags)
    if arm == :two_slice
        @test !any(outcome.syndrome)
        @test !any(outcome.recovery)
        if context.channel == "x_only"
            @test all(entry -> entry.logical_after == last(history).logical_after,history)
            @test Tuple(first(history).logical_after) == (outcome.logical_1,outcome.logical_2)
        end
    end
end

function audit_record(context,record)
    before = (copy(record.seed),copy(record.x),copy(record.z))
    outcomes = AF.three_outcomes(context,record;record_history=true)
    @test keys(outcomes) == (:final_only,:one_step,:two_slice)
    @test (record.seed,record.x,record.z) == before
    for (arm,outcome) in pairs(outcomes)
        audit_outcome(context,record,arm,outcome)
    end
    outcomes
end

if isdefined(AF,:build_context)
@testset "zero noise and final checkpoint conventions" begin
    for channel in ("x_only","z_only"), T in (1,2)
        context = AF.build_context(AF.AxisCode(T,2,2),0.0,channel)
        record = AF.sample_record(MersenneTwister(11),context)
        @test !any(record.seed) && !any(record.x) && !any(record.z)
        outcomes = audit_record(context,record)
        @test all(outcome -> !outcome.failure,values(outcomes))
        @test length(outcomes.final_only.history) == (channel == "z_only" ? 1 : 0)
        @test length(outcomes.one_step.history) == T
        @test length(outcomes.two_slice.history) == T+(channel == "z_only")
    end
end

@testset "every active postgate single fault and noisy seed" begin
    for channel in ("x_only","z_only")
        context = AF.build_context(AF.AxisCode(2,2,2),0.05,channel)
        locations = 0
        for (index,operation) in enumerate(context.encoder.operations), slot in eachindex(operation.qubits)
            record = cleanrecord(context)
            getproperty(record,channel == "x_only" ? :x : :z)[slot,index] = true
            audit_record(context,record)
            locations += 1
        end
        @test locations == 13context.code.T*context.code.Ly*context.code.Lz
        if channel == "z_only"
            for q in eachindex(cleanrecord(context).seed)
                record = cleanrecord(context)
                record.seed[q] = true
                audit_record(context,record)
            end
        end
    end
end

@testset "fixed random L3 records and immutable paired replay" begin
    for channel in ("x_only","z_only")
        context = AF.build_context(3,0.12,channel)
        rng = MersenneTwister(719)
        for _ in 1:12
            audit_record(context,AF.sample_record(rng,context))
        end
        a = AF.sample_record(MersenneTwister(13),context)
        b = AF.sample_record(MersenneTwister(13),context)
        @test (a.seed,a.x,a.z) == (b.seed,b.x,b.z)
        @test context.p == 0.12 && context.channel == channel
        if channel == "x_only"
            @test !any(a.seed) && !any(a.z)
        else
            @test !any(a.x)
        end
    end
end

@testset "seed repair preserves unobserved logical sectors" begin
    context = AF.build_context(3,0.05,"z_only")
    baseline = audit_record(context,cleanrecord(context))
    for logical in 1:2
        record = cleanrecord(context)
        record.seed[AF.logical_z(AF.AxisCode(0,3,3))[logical]] .= true
        outcomes = audit_record(context,record)
        for arm in keys(outcomes)
            outcome = outcomes[arm]
            @test outcome.logical_1 == (logical == 1)
            @test outcome.logical_2 == (logical == 2)
            @test !any(first(outcome.history).recovery)
            @test [entry.pre_syndrome for entry in outcome.history] ==
                [entry.pre_syndrome for entry in baseline[arm].history]
            @test [entry.recovery for entry in outcome.history] ==
                [entry.recovery for entry in baseline[arm].history]
        end
    end
end

@testset "one-edge endpoint route and syndrome-only decoder interface" begin
    context = AF.build_context(3,0.05,"x_only")
    kernel = context.slice_kernel
    singleton = falses(size(kernel.matrix,2)); singleton[1] = true
    syndrome = AF._axis_matrix_syndrome(kernel.matrix,singleton)
    decoded = AF._axis_endpoint_route(context,syndrome)
    @test length(only(decoded.paths)) == 1
    @test count(decoded.route) == 1
    @test AF._axis_matrix_syndrome(kernel.matrix,decoded.route) == syndrome
    @test AF._axis_decode(kernel,syndrome) == AF._axis_decode(kernel,copy(syndrome))
    @test_throws DimensionMismatch AF._axis_decode(kernel,falses(length(syndrome)+1))
    # A full-prefix X logical changes no plaquette measurement. Recovery must
    # therefore be identical while the final logical flag remains different.
    prefix = AF.AxisCode(1,3,3)
    firstwindow = first(context.windows)
    x0 = falses(AF.nqubits(context.code)); x1 = copy(x0)
    x1[AF.logical_x(prefix)[1]] .= true
    s0,s1 = parity(firstwindow.checks,x0),parity(firstwindow.checks,x1)
    @test s0 == s1
    @test AF._axis_decode(firstwindow.kernel,s0) == AF._axis_decode(firstwindow.kernel,s1)
    @test parity(AF.logical_z(context.code),x0) != parity(AF.logical_z(context.code),x1)
end

@testset "literal copy-stage feedforward and odd-parity retiring plane" begin
    context = AF.build_context(AF.AxisCode(2,2,3),0.05,"x_only")
    record = cleanrecord(context)
    for r in 1:context.code.T, (old,fresh) in zip(
            AF.plane_edges(context.code,r-1),AF.plane_edges(context.code,r))
        x,z = falses(AF.nqubits(context.code)),falses(AF.nqubits(context.code))
        x[old] = true
        firstop = r == 1 ? 1 : context.encoder.checkpoints[r-1]+1
        handsegment!(x,z,context,record,firstop,context.encoder.checkpoints[r])
        @test findall(x) == sort([old,fresh])
        @test !any(z)
    end
    context = AF.build_context(AF.AxisCode(2,2,3),0.05,"z_only")
    frame = AF.TN.PauliFrame(falses(AF.nqubits(context.code)),falses(AF.nqubits(context.code)))
    frame.z[AF.edge(context.code,:x,1,0,0)] = true
    entry = AF._axis_apply_z_push!(frame,context,1;record_history=true)
    @test count(entry.pre_syndrome) == 1
    @test count(entry.recovery) == 1
    @test !any(entry.prefix_syndrome)
    # Every allowed old x link retains its lower physical endpoint row.
    window = context.windows[2]
    for q in AF.slice_qubits(context.code,1)[1:context.code.Ly*context.code.Lz]
        @test count(check -> q in check,window.checks) == 2
    end
end

@testset "record validation and decoder reuse" begin
    for channel in ("x_only","z_only")
        context = AF.build_context(2,0.05,channel)
        record = cleanrecord(context)
        @test_throws DimensionMismatch AF.three_outcomes(context,AF.AxisRecord(falses(1),record.x,record.z))
        @test_throws DimensionMismatch AF.three_outcomes(context,AF.AxisRecord(record.seed,falses(2,1),record.z))
        bad = cleanrecord(context)
        getproperty(bad,channel == "x_only" ? :z : :x)[1,1] = true
        @test_throws ArgumentError AF.three_outcomes(context,bad)
        h = findfirst(operation -> operation.gate == :H,context.encoder.operations)
        bad = cleanrecord(context)
        getproperty(bad,channel == "x_only" ? :x : :z)[2,h] = true
        @test_throws ArgumentError AF.three_outcomes(context,bad)
        refreshed = AF.build_context(2,0.07,channel;previous=context)
        @test refreshed.p == 0.07
        @test refreshed.encoder === context.encoder
        @test all(window -> window.kernel === refreshed.windows[2].kernel,refreshed.windows[2:end])
        if channel == "z_only"
            @test refreshed.final === context.final
            @test refreshed.seed_kernel === context.seed_kernel
        else
            @test refreshed.final !== context.final
        end
    end
    @test_throws ArgumentError AF.build_context(2,-0.1,"z_only")
    @test_throws ArgumentError AF.build_context(2,0.1,"Y")
    @test_throws ArgumentError AF.build_context(AF.AxisCode(0,2,2),0.1,"x_only")
end
end
