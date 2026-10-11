using Test
using Random

@testset "3D scan can select only the Z channel" begin
    scan = compare_toric3d_feedback(MersenneTwister(42);
        sizes=[(2, 2, 2)], error_rates=[0.0], shots=1,
        seed=42, x_backend=:exact, channels=[:z_only])
    @test length(scan.points) == 1
    @test only(scan.points).channel === :z_only
    @test_throws ArgumentError compare_toric3d_feedback(MersenneTwister(42);
        sizes=[(2, 2, 2)], error_rates=[0.0], shots=1,
        seed=42, x_backend=:exact, channels=Symbol[])
end
