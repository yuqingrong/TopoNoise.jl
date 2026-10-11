using Test
using Random

@testset "3D paired comparison API" begin
    @test isdefined(TopoNoise, :compare_toric3d_feedback)
    @test isdefined(TopoNoise, :save_toric3d_comparison)
end

if isdefined(TopoNoise, :compare_toric3d_feedback)
    @testset "Paired estimates preserve discordance and finite uncertainty" begin
        m = TopoNoise._toric3d_metric(4, 2, 2, 1)
        @test m.final_only_failures == 2
        @test m.layer_feedback_failures == 2
        @test m.final_only_only == 1
        @test m.layer_feedback_only == 1
        @test m.difference == 0
        @test m.difference_standard_error ≈ sqrt(1 / 6)
        @test m.difference_ci[1] < 0 < m.difference_ci[2]
        z = TopoNoise._toric3d_metric(1000, 0, 0, 0)
        @test z.final_only_ci[1] == 0
        @test 0 < z.final_only_ci[2] < .004
        @test z.difference_ci[1] < 0 < z.difference_ci[2]
        @test toric3d_z_lower_bound(.01, 2, 2) ≈ (1 - .98^4) / 2
    end

    @testset "Seeded scan and checkpoint recovery" begin
        args = (sizes=[(2, 2, 2)], error_rates=[0., .03], shots=12,
                seed=77, x_backend=:exact)
        a = compare_toric3d_feedback(MersenneTwister(77); args...)
        b = compare_toric3d_feedback(MersenneTwister(77); args...)
        @test length(a.points) == 4
        @test a.dimensions == [(2, 2, 2)]
        for (x, y) in zip(a.points, b.points)
            @test x.metrics == y.metrics
            @test x.p == y.p
            @test x.channel == y.channel
            @test x.shots == 12
            for metric in values(x.metrics)
                @test 0 <= metric.final_only_failures <= 12
                @test 0 <= metric.layer_feedback_failures <= 12
            end
            @test x.metrics.state_failure == x.metrics.logical_x
            if x.p == 0
                @test all(m.final_only_failures == m.layer_feedback_failures == 0
                          for m in values(x.metrics))
            end
            inactive = x.channel === :x_only ? x.metrics.logical_z : x.metrics.logical_x
            @test inactive.final_only_failures == inactive.layer_feedback_failures == 0
        end
        @test_throws ArgumentError compare_toric3d_feedback(MersenneTwister(1); sizes=[], shots=1)
        @test_throws ArgumentError compare_toric3d_feedback(MersenneTwister(1); shots=0)
        @test_throws ArgumentError compare_toric3d_feedback(MersenneTwister(1); error_rates=[.1, .01])
        @test_throws ArgumentError compare_toric3d_feedback(MersenneTwister(1); sizes=[2,2])

        mktempdir() do dir
            c = compare_toric3d_feedback(MersenneTwister(77); args..., output_dir=dir)
            raw = joinpath(dir, "toric3d-feedback-raw.csv")
            @test isfile(raw)
            lines = readlines(raw)
            @test length(lines) == 1 + 4 * length(c.points)
            @test occursin("difference_feedback_minus_final", first(lines))
            @test occursin("pymatching_version", first(lines))
            @test isfile(joinpath(dir, "toric3d-feedback-metadata.toml"))
            paths = save_toric3d_comparison(c, dir)
            for key in (:csv, :metadata, :svg, :pdf, :png)
                @test isfile(getproperty(paths, key))
                @test filesize(getproperty(paths, key)) > 0
            end
        end
    end

    @testset "Layer feedback equals propagated correction frame" begin
        encoder = toric3d_encoder(ToricCode3D(2,2,2))
        decoders = [build_toric3d_decoder(code; p=.03, x_backend=:exact)
                    for code in encoder.prefixes]
        n = data_qubit_count(encoder.code)
        for seed in 1:10
            faults = sample_toric3d_faults(MersenneTwister(seed), encoder; p_x=.03,p_z=.04)
            actual = replay_toric3d_feedback(encoder,faults,decoders; record_history=true)
            uncorrected = PauliFrame(falses(n),falses(n))
            deferred = PauliFrame(falses(n),falses(n))
            first_gate = 1
            for j in eachindex(encoder.checkpoints)
                last_gate = encoder.checkpoints[j]
                propagate_toric3d_segment!(uncorrected,encoder,faults,first_gate,last_gate)
                negative = Ref(false)
                for operation in encoder.operations[first_gate:last_gate]
                    TopoNoise._conjugate_pauli!(deferred.x,deferred.z,negative,operation)
                end
                mapping = encoder.prefix_qubits[j]
                effective = PauliFrame((uncorrected.x .⊻ deferred.x)[mapping],
                                       (uncorrected.z .⊻ deferred.z)[mapping])
                s = toric3d_syndrome(encoder.prefixes[j],effective)
                @test s == actual.history[j].syndrome
                recovery = decode_toric3d(decoders[j],s.stars,s.plaquettes)
                deferred.x[mapping] .⊻= recovery.x
                deferred.z[mapping] .⊻= recovery.z
                first_gate = last_gate + 1
            end
            @test actual.frame == PauliFrame(uncorrected.x .⊻ deferred.x,
                                             uncorrected.z .⊻ deferred.z)
        end
    end

    @testset "Production BP+OSD scan is reproducible and records its version" begin
        kwargs = (sizes=[(1,2,2)], error_rates=[.03], shots=16, seed=918)
        a = compare_toric3d_feedback(MersenneTwister(918);kwargs...)
        b = compare_toric3d_feedback(MersenneTwister(918);kwargs...)
        @test [p.metrics for p in a.points] == [p.metrics for p in b.points]
        mktempdir() do dir
            paths = TopoNoise._write_toric3d_checkpoint(a,dir)
            lines = readlines(paths.csv)
            header = split(lines[1],',')
            column = findfirst(==("ldpc_version"),header)
            @test column !== nothing
            @test all(strip(split(line,',')[column], '"') == first(a.points).decoder_metadata.x_version
                      for line in lines[2:end])
        end
    end

    @testset "First-layer Z pair is indistinguishable even with feedback" begin
        encoder = toric3d_encoder(ToricCode3D(2,2,2))
        decoders = [build_toric3d_decoder(code; p=.03, x_backend=:exact)
                    for code in encoder.prefixes]
        for seed in 0:5
            base = sample_toric3d_faults(MersenneTwister(seed),encoder;
                p_x=seed == 0 ? 0.0 : .02,p_z=seed == 0 ? 0.0 : .03)
            a = with_toric3d_fault(encoder,base,2,data_qubit_index(encoder.code,:x,0,0,0),:Z)
            b = with_toric3d_fault(encoder,base,6,data_qubit_index(encoder.code,:x,1,0,0),:Z)
            for protocol in (:final_only,:layer_feedback)
                ra = replay_toric3d_feedback(encoder,a,decoders;protocol,record_history=true)
                rb = replay_toric3d_feedback(encoder,b,decoders;protocol,record_history=true)
                @test [h.syndrome for h in ra.history] == [h.syndrome for h in rb.history]
                @test ra.logicals.logical_z != rb.logicals.logical_z
                @test ra.logicals.logical_x == rb.logicals.logical_x
            end
        end
    end
end
