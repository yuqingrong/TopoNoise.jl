using Test, Random, TopoNoise
import TOML, SHA

const FULL_PREFIX_RUNNER_PATH=joinpath(@__DIR__,"run_study.jl")
@testset "H+6 paired full-prefix X/Z runner exists" begin
    @test isfile(FULL_PREFIX_RUNNER_PATH)
end

if isfile(FULL_PREFIX_RUNNER_PATH)
include(FULL_PREFIX_RUNNER_PATH)
@testset "Full-prefix runner protocol binds every-transfer correction" begin
    @test H6FullPrefixStudy.ARMS == ("final_only","full_prefix")
    for channel in ("x_only","z_only")
        protocol=H6FullPrefixStudy.protocol_metadata(channel)
        @test protocol["feedback_timing"] == "after transfers 1:T inclusive; no separate r=0 correction"
        @test occursin("Q0..Qr",protocol["full_prefix"])
        @test occursin("family.prefixes[end]",protocol["shared_final_decoder"])
    end
end

function full_prefix_args(directory;channel="x_only",rates="0,.05",shots="4",seed="1234",extra=String[])
    return vcat(["--channel",channel,"--output",directory,"--sizes","2","--rates",rates,
        "--shots",shots,"--seed",seed,"--batch-size","2"],extra)
end

function full_prefix_raw_counts(path)
    rows=split.(readlines(path),',')
    # Column 11 is elapsed wall time; every other field must match bit-for-bit.
    return [row[[1:10;12:length(row)]] for row in rows]
end

@testset "Full-prefix CLI, point seeds, and four-bin encoding" begin
    options=H6FullPrefixStudy.parse_arguments(full_prefix_args("out"))
    @test options.config["channel"] == "x_only"
    @test options.config["sizes"] == [2]
    @test options.config["rates"] == [0.0,0.05]
    @test !options.resume
    @test options.max_batches == typemax(Int)
    @test H6FullPrefixStudy.failure_mask((false,false)) == 0
    @test H6FullPrefixStudy.failure_mask((true,false)) == 1
    @test H6FullPrefixStudy.failure_mask((false,true)) == 2
    @test H6FullPrefixStudy.failure_mask((true,true)) == 3
    @test H6FullPrefixStudy.point_seeds(1234,2) == UInt64[0xe4974282f47c0920,0xfaada3ebf6911431]
    @test length(unique(H6FullPrefixStudy.point_seeds(1234,20))) == 20
    for extra in (["--channel","y_only"],["--sizes","2,2"],["--sizes","1"],["--rates",".05,.01"],
            ["--rates",".05,.05"],["--rates","-.01"],["--rates",".5"],
            ["--shots","0"],["--seed","-1"],["--batch-size","0"],
            ["--max-batches","0"],["--unknown","1"])
        args=full_prefix_args("out")
        index=findfirst(==(extra[1]),args)
        index === nothing ? append!(args,extra) : (args[index+1]=extra[2])
        @test_throws ArgumentError H6FullPrefixStudy.parse_arguments(args)
    end
    @test_throws ArgumentError H6FullPrefixStudy.parse_arguments(["--channel"])
end

@testset "Full-prefix exact next-record resume and failure diagnostics" begin
    encoder=toric3d_encoded_encoder(Toric3DEncodedSlice(2,2,2))
    for channel in ("x_only","z_only")
        rng=Random.Xoshiro(UInt64(0x1234))
        sample(rng)=channel == "x_only" ? sample_toric3d_faults(rng,encoder;p_x=.05,p_z=0.) :
            sample_toric3d_diagonal_z_faults(rng,encoder;p=.05)
        records=[sample(rng) for _ in 1:4]
        for completed in 0:3
            restored=H6FullPrefixStudy.resume_rng(UInt64(0x1234),completed,encoder,.05,channel)
            @test sample(restored) == records[completed+1]
        end
        mktempdir() do directory
            H6FullPrefixStudy.dump_failure(directory,(L=2,p=.05,channel=channel),UInt64(0x1234),3,records[3],ErrorException("invariant fixture"))
            dump=TOML.parsefile(joinpath(directory,"failed-record.toml"))
            @test dump["channel"] == channel
            @test dump["point_seed_hex"] == "0000000000001234"
            @test dump["shot_index"] == 3
            @test dump["dimensions"] == [2,2,2]
            transfer=channel == "x_only" ? records[3] : records[3].transfer
            @test dump["transfer_shape"] == collect(size(transfer.z))
            @test dump["transfer_x_true_indices"] == findall(vec(transfer.x))
            @test dump["transfer_z_true_indices"] == findall(vec(transfer.z))
            if channel == "z_only"
                @test dump["seed_true_indices"] == findall(records[3].seed)
            end
        end
    end
end

@testset "Full-prefix X/Z paired scans, channel caches, exact partial resume" begin
    for channel in ("x_only","z_only")
        mktempdir() do directory
            direct=joinpath(directory,"direct"); resumed=joinpath(directory,"resumed")
            partial=H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel,extra=["--max-batches","3"]));io=devnull)
            @test partial["completed_shots"] == 6
            @test !partial["complete"]
            state=TOML.parsefile(joinpath(resumed,"checkpoint.toml"))
            @test [s["completed_shots"] for s in state["states"]] == [4,2]
            H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel,extra=["--resume"]));io=devnull)
            full=H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(full_prefix_args(direct;channel));io=devnull)
            @test full["complete"]
            @test full["completed_shots"] == 8
            @test length(readlines(joinpath(direct,"raw.csv"))) == 3
            # X priors depend on p; Z matching weights do not.
            @test full["decoder_build_count_by_size"]["2"] == (channel == "x_only" ? 2 : 1)
            @test full["decoder_build_seconds_by_size"]["2"] >= 0
            @test full_prefix_raw_counts(joinpath(direct,"raw.csv")) == full_prefix_raw_counts(joinpath(resumed,"raw.csv"))
            state=TOML.parsefile(joinpath(resumed,"checkpoint.toml"))
            @test state["states"][1]["joint_counts"] == [4,0,0,0]
            @test state["states"][1]["failures"] == [0,0]
            for point in state["states"]
                @test sum(point["joint_counts"]) == point["completed_shots"]
                @test point["failures"][1] == point["joint_counts"][2]+point["joint_counts"][4]
                @test point["failures"][2] == point["joint_counts"][3]+point["joint_counts"][4]
            end
            raw_before=read(joinpath(resumed,"raw.csv"))
            checkpoint_before=read(joinpath(resumed,"checkpoint.toml"))
            complete=H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel,extra=["--resume"]));io=devnull)
            @test complete["complete"]
            @test read(joinpath(resumed,"raw.csv")) == raw_before
            @test read(joinpath(resumed,"checkpoint.toml")) == checkpoint_before
            @test complete["raw_sha256"] == bytes2hex(SHA.sha256(raw_before))
            @test_throws ArgumentError H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(full_prefix_args(resumed;channel));io=devnull)
            @test_throws ArgumentError H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel,shots="5",extra=["--resume"]));io=devnull)
            other_channel=channel == "x_only" ? "z_only" : "x_only"
            @test_throws ArgumentError H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel=other_channel,extra=["--resume"]));io=devnull)
            state["states"][2]["joint_counts"][1] += 1
            open(joinpath(resumed,"checkpoint.toml"),"w") do io
                TOML.print(io,state;sorted=true)
            end
            @test_throws ArgumentError H6FullPrefixStudy.run_study(H6FullPrefixStudy.parse_arguments(
                full_prefix_args(resumed;channel,extra=["--resume"]));io=devnull)
        end
    end
end
end
