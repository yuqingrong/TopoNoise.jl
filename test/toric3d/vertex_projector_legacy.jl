using Test, Random, TopoNoise

@testset "Legacy H+6 plaquette-transfer schedule is frozen" begin
    e = toric3d_slice_encoder(Toric3DSliceLattice(1,2,2))
    @test e.checkpoints == [28]
    @test sum(length(o.qubits) for o in e.operations) == 52
    @test [Tuple(o.qubits) for o in e.operations] == [
        (13,), (13,17), (13,21), (7,17), (2,17), (11,21), (1,21),
        (14,), (14,18), (14,22), (8,18), (1,18), (12,22), (2,22),
        (15,), (15,19), (15,23), (5,19), (4,19), (9,23), (3,23),
        (16,), (16,20), (16,24), (6,20), (3,20), (10,24), (4,24)]
    @test [o.gate for o in e.operations] == repeat([:H,fill(:CNOT,6)...],4)
end

@testset "Legacy H+6 replay and artifact formats remain reproducible" begin
    c=Toric3DEncodedSlice(1,2,2)
    e=toric3d_encoded_encoder(c)
    r=sample_toric3d_diagonal_z_faults(MersenneTwister(17),e;p=.23)
    d=build_toric3d_z_feedback_decoders(c)
    a=toric3d_z_strength_triplet(e,r,d)
    @test r == sample_toric3d_diagonal_z_faults(MersenneTwister(17),e;p=.23)
    @test [(a[k].logical_1,a[k].logical_2,a[k].failure) for k in keys(a)] ==
        [(false,true,true),(true,true,true),(true,true,true)]
    # The existing independent literal header/CSV fixture is the schema oracle.
    include("feedback_strength_output.jl")
end
