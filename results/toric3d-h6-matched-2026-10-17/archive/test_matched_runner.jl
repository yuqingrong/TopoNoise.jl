using Test, TopoNoise, Random
import TOML

const RUNNER = joinpath(@__DIR__, "run_study.jl")

@testset "Matched runner implementation exists" begin
    @test isfile(RUNNER)
end

include(RUNNER)
const Study = H6MatchedStudy

function compare_outcome(actual, expected)
    @test actual.frame.x == expected.frame.x
    @test actual.frame.z == expected.frame.z
    for key in (:syndrome, :recovery, :residual, :logical_1, :logical_2, :failure)
        @test getproperty(actual, key) == getproperty(expected, key)
    end
end

function options(directory, channel; resume=false, max_batches=typemax(Int))
    args = ["--channel", channel, "--output", directory, "--sizes", "3",
        "--rates", "0,0.15", "--shots", "5", "--seed", "20261017",
        "--batch-size", "2", "--max-batches", string(max_batches)]
    resume && push!(args, "--resume")
    return Study.parse_arguments(args)
end

@testset "Fixed fault records reach each existing protocol unchanged" begin
    for channel in ("x_only", "z_only")
        context = Study.build_context(3, 0.15, channel)
        rng = Random.Xoshiro(0x91872)
        for p in (0.0, 0.15, 0.2)
            raw = Study.sample_faults(rng, context.encoder, p, channel)
            saved = deepcopy(raw)
            actual = Study.three_outcomes(context, raw, channel; record_history=true)
            @test raw == saved
            @test length(actual) == 3
            @test all(result.failure == (result.logical_1 | result.logical_2) for result in actual)
            if channel == "x_only"
                expected = [replay_toric3d_x_feedback(context.encoder, raw, context.slice, context.two.final;
                        feedback=false, routing=:one_step, record_history=true),
                    replay_toric3d_x_feedback(context.encoder, raw, context.slice, context.two.final;
                        feedback=true, routing=:one_step, record_history=true),
                    replay_toric3d_h6_two_slice_x(context.encoder, raw, context.two;
                        feedback=true, record_history=true)]
                for arm in 1:3
                    compare_outcome(actual[arm], expected[arm])
                end
            else
                corrected, syndrome, recovery = Study.common_seed_record(raw, context.family)
                corrected_saved = deepcopy(corrected)
                @test corrected.transfer == raw.transfer
                @test corrected.seed == raw.seed .⊻ recovery
                @test !any(isodd.(context.family.seed.matrix * UInt8.(corrected.seed)))
                original_one = replay_toric3d_h6_z_one_step(context.encoder, raw, context.one;
                    record_history=true)
                compare_outcome(actual[2], original_one)
                @test syndrome == first(original_one.history).syndrome
                @test recovery == first(original_one.history).recovery
                compare_outcome(actual[1], replay_toric3d_h6_z_one_step(context.encoder, raw,
                    context.one; intermediate_feedback=false))
                expected = [replay_toric3d_h6_two_slice_z(context.encoder, corrected, context.family;
                        feedback=false, record_history=true),
                    replay_toric3d_h6_z_one_step(context.encoder, corrected, context.one;
                        record_history=true),
                    replay_toric3d_h6_two_slice_z(context.encoder, corrected, context.family;
                        feedback=true, record_history=true)]
                @test corrected == corrected_saved
                @test raw == saved
                @test !any(first(actual[2].history).syndrome)
                @test !any(first(actual[2].history).recovery)
                for arm in 1:3
                    compare_outcome(actual[arm], expected[arm])
                end
            end
            p == 0 && @test all(!result.failure for result in actual)
        end
    end
end

@testset "Common seed repair retains undetectable seed logical errors" begin
    context = Study.build_context(3, 0.0, "z_only")
    clean = Study.sample_faults(Random.Xoshiro(9), context.encoder, 0.0, "z_only")
    for support in toric3d_encoded_logical_z(context.code)
        seed = copy(clean.seed)
        seed[support] .= true
        raw = Toric3DDiagonalZFaultRecord(seed, clean.transfer)
        corrected, syndrome, recovery = Study.common_seed_record(raw, context.family)
        @test !any(syndrome)
        @test !any(recovery)
        @test corrected.seed == seed
        @test all(result.failure for result in Study.three_outcomes(context, raw, "z_only"))
    end
end

@testset "Three marginals are exactly the eight shared-record bins" begin
    config = options("/private/tmp/unused-matched-test", "x_only").config
    config["shots"] = 8
    checkpoint = Study.fresh_checkpoint(config, Dict(), Dict())
    state = first(checkpoint["states"])
    for mask in 0:7
        outcomes = [(logical_1=(mask & (1 << (arm-1))) != 0, logical_2=false,
            failure=(mask & (1 << (arm-1))) != 0) for arm in 1:3]
        Study.accumulate!(state, outcomes)
    end
    @test state["completed_shots"] == 8
    @test state["joint_counts"] == ones(Int, 8)
    @test state["failures"] == [4, 4, 4]
    @test state["logical_1_counts"] == [4, 4, 4]
    @test state["logical_2_counts"] == [0, 0, 0]
    @test Study.validate_states(checkpoint) === checkpoint
    for mutate in (
            c -> c["states"][1]["joint_counts"][1] += 1,
            c -> c["states"][1]["failures"][2] -= 1,
            c -> c["states"][1]["logical_2_counts"][3] = 8,
            c -> c["point_seeds_hex"][1] = "0000000000000000",
            c -> c["states"][2]["completed_shots"] = true,
            c -> c["states"][2]["elapsed_seconds"] = NaN,
            c -> c["generation"] = -1,
            c -> c["updated_at_unix"] = NaN)
        damaged = deepcopy(checkpoint)
        mutate(damaged)
        @test_throws ArgumentError Study.validate_states(damaged)
    end
end

@testset "Checkpoint resume consumes exactly the original point streams" begin
    for channel in ("x_only", "z_only")
        mktempdir() do directory
            full = joinpath(directory, "full")
            resumed = joinpath(directory, "resumed")
            Study.run_study(options(full, channel); io=devnull)
            partial = Study.run_study(options(resumed, channel; max_batches=2); io=devnull)
            @test partial["completed_shots"] == 4
            @test !partial["complete"]
            # Mirrors can be stale or absent: the authoritative checkpoint repairs them.
            write(joinpath(resumed, "raw.csv"), "stale mirror")
            rm(joinpath(resumed, "metadata.toml"))
            done = Study.run_study(options(resumed, channel; resume=true); io=devnull)
            @test done["complete"]
            uninterrupted = TOML.parsefile(joinpath(full, "checkpoint.toml"))
            restarted = TOML.parsefile(joinpath(resumed, "checkpoint.toml"))
            @test uninterrupted["point_seeds_hex"] == restarted["point_seeds_hex"]
            for (a, b) in zip(uninterrupted["states"], restarted["states"])
                for key in ("completed_shots", "logical_1_counts", "logical_2_counts", "failures", "joint_counts")
                    @test a[key] == b[key]
                end
            end
            @test length(readlines(joinpath(resumed, "raw.csv"))) == 3
            @test occursin("one_step_failures", first(readlines(joinpath(resumed, "raw.csv"))))
            @test occursin("joint_07", first(readlines(joinpath(resumed, "raw.csv"))))
            @test isfile(joinpath(resumed, "previous-checkpoint.toml"))
            @test !isdir(joinpath(resumed, "checkpoints"))
            # Checksum and provenance mismatches must reject, never use the valid backup.
            original = read(joinpath(resumed, "checkpoint.toml"), String)
            altered = deepcopy(restarted)
            altered["states"][1]["failures"][1] += 1
            write(joinpath(resumed, "checkpoint.toml"), Study.toml_text(altered))
            @test_throws ArgumentError Study.run_study(options(resumed, channel; resume=true); io=devnull)
            write(joinpath(resumed, "checkpoint.toml"), original)
            config = options(resumed, channel).config
            bad_config = deepcopy(config); bad_config["shots"] += 1
            @test_throws ArgumentError Study.load_checkpoint(resumed, bad_config,
                restarted["runtime"], restarted["source_hashes"])
            @test_throws ArgumentError Study.load_checkpoint(resumed, config,
                merge(restarted["runtime"], Dict("julia_version"=>"changed")), restarted["source_hashes"])
            @test_throws ArgumentError Study.load_checkpoint(resumed, config,
                restarted["runtime"], merge(restarted["source_hashes"], Dict("run_study.jl"=>"changed")))
            backup_path = joinpath(resumed, "previous-checkpoint.toml")
            backup_text = read(backup_path, String)
            corrupted_backup = TOML.parse(backup_text)
            corrupted_backup["generation"] += 1
            write(backup_path, Study.toml_text(corrupted_backup))
            @test_throws ArgumentError Study.load_checkpoint(resumed, config,
                restarted["runtime"], restarted["source_hashes"])
            write(backup_path, backup_text)
            # A malformed authoritative file is rejected even with a valid backup.
            write(joinpath(resumed, "checkpoint.toml"), "[unterminated")
            @test_throws ArgumentError Study.run_study(options(resumed, channel; resume=true); io=devnull)
        end
    end
end

@testset "Source fingerprints do not depend on checkout location" begin
    mktempdir() do directory
        for name in ("a", "b")
            root = joinpath(directory, name, "runtime")
            mkpath(joinpath(root, "src", "nested"))
            write(joinpath(root, "Project.toml"), "name = \"Fixture\"\n")
            write(joinpath(root, "Manifest.toml"), "julia_version = \"1.10.12\"\n")
            write(joinpath(root, "src", "nested", "a.jl"), "const A=1\n")
            write(joinpath(root, "src", "data.txt"), "included non-Julia source asset\n")
            write(joinpath(directory, name, "run_study.jl"), "fixture runner\n")
        end
        a = joinpath(directory, "a"); b = joinpath(directory, "b")
        hashes = Study.source_hashes(joinpath(a, "runtime"), joinpath(a, "run_study.jl"))
        @test hashes == Study.source_hashes(joinpath(b, "runtime"), joinpath(b, "run_study.jl"))
        write(joinpath(b, "runtime", "src", "data.txt"), "changed asset\n")
        @test hashes != Study.source_hashes(joinpath(b, "runtime"), joinpath(b, "run_study.jl"))
    end
end
