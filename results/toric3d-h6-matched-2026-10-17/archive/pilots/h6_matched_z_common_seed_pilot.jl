using TopoNoise, Random, LinearAlgebra
import PythonCall

function common_seed_record(raw, family)
    @assert family.seed.window.qubits == toric3d_slice_qubits(family.code.lattice,0)
    syndrome = BitVector(isodd.(family.seed.matrix * UInt8.(raw.seed)))
    recovery = decode_toric3d_z_feedback(family.seed, syndrome)
    seed = raw.seed .⊻ recovery
    @assert !any(isodd.(family.seed.matrix * UInt8.(seed)))
    return Toric3DDiagonalZFaultRecord(seed,raw.transfer), syndrome,recovery
end

function outcomes(encoder,raw,family,one_step; history=false)
    corrected,s,c = common_seed_record(raw,family)
    baseline = replay_toric3d_h6_two_slice_z(encoder,corrected,family;feedback=false,record_history=history)
    one = replay_toric3d_h6_z_one_step(encoder,corrected,one_step;record_history=history)
    two = replay_toric3d_h6_two_slice_z(encoder,corrected,family;record_history=history)
    return baseline,one,two,corrected,s,c
end

function compare_result(a,b)
    @assert a.frame.x == b.frame.x && a.frame.z == b.frame.z
    @assert all(getproperty(a,k)==getproperty(b,k) for k in (:syndrome,:recovery,:residual,:logical_1,:logical_2,:failure))
end

function checks()
    cases=0
    for dims in ((1,2,2),(3,2,3),(3,3,3))
        code=Toric3DEncodedSlice(dims...);encoder=toric3d_encoded_encoder(code)
        family=build_toric3d_z_feedback_decoders(code);one_step=build_toric3d_h6_z_one_step_decoder(family)
        rng=Random.Xoshiro(0x12345678)
        records=[sample_toric3d_diagonal_z_faults(rng,encoder;p=p) for p in (0.0,0.02,0.2)]
        clean=first(records)
        seed=copy(clean.seed);seed[1]=true
        push!(records,Toric3DDiagonalZFaultRecord(seed,clean.transfer))
        for support in toric3d_encoded_logical_z(code)
            seed=copy(clean.seed);seed[support].=true
            push!(records,Toric3DDiagonalZFaultRecord(seed,clean.transfer))
        end
        gate=length(encoder.operations)
        lastfault=with_toric3d_fault(encoder,clean.transfer,gate,last(encoder.operations[gate].qubits),:Z)
        push!(records,Toric3DDiagonalZFaultRecord(clean.seed,lastfault))
        for raw in records
            saved=deepcopy(raw)
            base,one,two,corrected,s,c=outcomes(encoder,raw,family,one_step;history=true)
            @assert raw==saved
            @assert corrected.transfer==raw.transfer
            @assert corrected.seed == raw.seed .⊻ c
            original_one=replay_toric3d_h6_z_one_step(encoder,raw,one_step;record_history=true)
            compare_result(one,original_one)
            compare_result(base,replay_toric3d_h6_z_one_step(encoder,raw,one_step;intermediate_feedback=false))
            @assert first(original_one.history).syndrome==s
            @assert first(original_one.history).recovery==c
            @assert first(one.history).kind===:seed && !any(first(one.history).syndrome) && !any(first(one.history).recovery)
            @assert length(one.history)==dims[1]
            @assert isempty(base.history)
            @assert length(two.history)==dims[1]
            @assert [h.checkpoint for h in two.history]==collect(1:dims[1])
            @assert all(h.kind===:window && !any(h.post_syndrome) for h in two.history)
            @assert !any(two.syndrome) && !any(two.recovery)
            for result in (base,one,two)
                @assert !any(isodd.(family.prefixes[end].matrix*UInt8.(result.residual)))
                @assert result.failure==(result.logical_1|result.logical_2)
            end
            cases+=1
        end
        for raw in records[5:6]
            @assert all(r.failure for r in outcomes(encoder,raw,family,one_step)[1:3])
        end
    end
    println("CHECKS_PASS cases=",cases);flush(stdout)
end

function main()
    @assert Threads.nthreads()==1 && BLAS.get_num_threads()==1
    println("RUNTIME Julia=",VERSION," Python=",PythonCall.pyconvert(String,PythonCall.pyimport("sys").executable));flush(stdout)
    checks()
    println("TIMING_HEADER,L,p,records,build_seconds,sample_s,seed_repair_s,final_only_s,one_step_s,two_slice_s,total_s")
    for L in (9,11,13,15)
        build_s=@elapsed begin
            code=Toric3DEncodedSlice(L,L,L);encoder=toric3d_encoded_encoder(code)
            family=build_toric3d_z_feedback_decoders(code)
            one_step=build_toric3d_h6_z_one_step_decoder(family)
        end
        rng=Random.Xoshiro(UInt64(2026101700+L))
        outcomes(encoder,sample_toric3d_diagonal_z_faults(rng,encoder;p=.02),family,one_step)
        for p in (.005,.02,.035)
            durations=zeros(5)
            elapsed=@elapsed for shot in 1:3
                durations[1]+=@elapsed raw=sample_toric3d_diagonal_z_faults(rng,encoder;p)
                saved=deepcopy(raw)
                durations[2]+=@elapsed (corrected,s,c)=common_seed_record(raw,family)
                durations[3]+=@elapsed a=replay_toric3d_h6_two_slice_z(encoder,corrected,family;feedback=false)
                @assert raw==saved
                durations[4]+=@elapsed b=replay_toric3d_h6_z_one_step(encoder,corrected,one_step)
                @assert raw==saved
                durations[5]+=@elapsed c=replay_toric3d_h6_two_slice_z(encoder,corrected,family)
                @assert raw==saved
            end
            println(join(["TIMING",L,p,3,build_s,durations...,elapsed],','));flush(stdout)
        end
        GC.gc()
    end
end
main()
