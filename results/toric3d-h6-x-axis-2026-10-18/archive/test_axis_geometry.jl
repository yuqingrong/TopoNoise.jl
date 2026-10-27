using Test

# Keep geometry tests independent of the concurrently developed feedback module.
module AxisGeometryUnderTest
using TopoNoise
const TN = TopoNoise
const source = joinpath(@__DIR__, "axis_geometry.jl")
isfile(source) && include(source)
end
const AG = AxisGeometryUnderTest

@testset "physical-x H+6 geometry API" begin
    @test isdefined(AG, :AxisCode)
    @test isdefined(AG, :axis_encoder)
    @test isdefined(AG, :verify_encoder)
end

function checkmatrix(checks, n)
    matrix = falses(length(checks), n)
    for (row, support) in enumerate(checks), qubit in support
        matrix[row, qubit] = !matrix[row, qubit]
    end
    matrix
end

# Independent elimination, so topology tests do not rely on the verifier's basis.
function rank2(input)
    matrix = copy(input)
    rank = 0
    for column in axes(matrix, 2)
        pivot = findfirst(row -> matrix[row, column], rank+1:size(matrix, 1))
        isnothing(pivot) && continue
        source = rank + pivot
        rank += 1
        matrix[rank, :], matrix[source, :] = copy(matrix[source, :]), copy(matrix[rank, :])
        for row in rank+1:size(matrix, 1)
            matrix[row, column] && (matrix[row, :] .⊻= matrix[rank, :])
        end
    end
    rank
end

@testset "square-seed slab topology and logical pairs" begin
    for (T, Ly, Lz) in ((0, 2, 2), (0, 3, 4), (1, 2, 2), (2, 2, 3),
                        (3, 3, 3), (2, 3, 4))
        code = AG.AxisCode(T, Ly, Lz)
        N = Ly * Lz
        n = AG.nqubits(code)
        @test n == (3T + 2)N
        @test vcat([AG.slice_qubits(code, r) for r in 0:T]...) == collect(1:n)
        for r in 0:T
            @test AG.plane_edges(code, r) == vcat(
                [AG.edge(code, :y, r, u, v) for u in 0:Ly-1 for v in 0:Lz-1],
                [AG.edge(code, :z, r, u, v) for u in 0:Ly-1 for v in 0:Lz-1])
            @test AG.edge(code, :y, r, -1, Lz) == AG.edge(code, :y, r, Ly-1, 0)
        end
        stars, plaquettes = AG.star_checks(code), AG.plaquette_checks(code)
        @test length(stars) == (T+1)N
        @test length(plaquettes) == (3T+1)N
        @test all(length(unique(check)) == length(check) for check in vcat(stars, plaquettes))
        @test all(length(check) == 4 for check in plaquettes)
        @test [length(stars[r*N+1]) for r in 0:T] ==
            [4 + (r > 0) + (r < T) for r in 0:T]
        A, B = checkmatrix(stars, n), checkmatrix(plaquettes, n)
        @test all(iseven, Int.(A) * Int.(B'))
        @test all(==(2), vec(sum(A; dims=1)))
        @test rank2(A) == (T+1)N - 1
        @test rank2(B) == (2T+1)N - 1
        @test n - rank2(A) - rank2(B) == 2
        X = checkmatrix(AG.logical_x(code), n)
        Z = checkmatrix(AG.logical_z(code), n)
        @test all(iseven, Int.(X) * Int.(B'))
        @test all(iseven, Int.(Z) * Int.(A'))
        @test mod.(Int.(X) * Int.(Z'), 2) == [1 0; 0 1]
        @test length.(AG.logical_z(code)) == [Ly, Lz]
        @test length.(AG.logical_x(code)) == [(T+1)Lz, (T+1)Ly]
    end
    @test_throws ArgumentError AG.AxisCode(-1, 3, 3)
    @test_throws ArgumentError AG.AxisCode(2, 1, 3)
    @test_throws ArgumentError AG.AxisCode(2, 3, 1)
    code = AG.AxisCode(2, 2, 3)
    @test_throws ArgumentError AG.edge(code, :x, 0, 0, 0)
    @test_throws ArgumentError AG.edge(code, :y, 3, 0, 0)
    @test_throws ArgumentError AG.edge(code, :a, 1, 0, 0)
    @test_throws ArgumentError AG.slice_qubits(code, -1)
    # The two period-two edges share endpoints but are distinct physical qubits.
    @test AG.edge(code, :y, 0, 0, 0) != AG.edge(code, :y, 0, 1, 0)
end

@testset "H+6 schedule, noise locations, and causal targets" begin
    for (T, Ly, Lz) in ((0, 2, 3), (1, 2, 2), (2, 2, 3), (3, 3, 3))
        code = AG.AxisCode(T, Ly, Lz)
        encoder = AG.axis_encoder(code)
        N = Ly * Lz
        @test encoder.code == code
        @test encoder.checkpoints == [7N*r for r in 1:T]
        @test length(encoder.operations) == 7N*T
        for r in 1:T
            ops = encoder.operations[7N*(r-1)+1:7N*r]
            fresh = Set(AG.slice_qubits(code, r))
            @test count(op -> op.gate === :H, ops) == N
            @test count(op -> op.gate === :CNOT, ops) == 6N
            @test sum(length(op.qubits) for op in ops) == 13N
            oldplane, newplane = AG.plane_edges(code, r-1), AG.plane_edges(code, r)
            @test [op.qubits for op in ops[1:2N]] == [[a,b] for (a,b) in zip(oldplane,newplane)]
            @test all(op.gate === :CNOT for op in ops[1:2N])
            for (i, (u,v)) in enumerate((u,v) for u in 0:Ly-1 for v in 0:Lz-1)
                vertexops = ops[2N+5(i-1)+1:2N+5i]
                link = AG.edge(code, :x, r, u, v)
                @test vertexops[1].gate === :H
                @test vertexops[1].qubits == [link]
                @test [op.qubits for op in vertexops[2:5]] ==
                    [[link, AG.edge(code, d, r, a, b)]
                     for (d,a,b) in ((:y,u,v),(:y,u-1,v),(:z,u,v),(:z,u,v-1))]
            end
            @test all(last(op.qubits) in fresh for op in ops)
            @test all(all(q <= AG.nqubits(AG.AxisCode(r,Ly,Lz)) for q in op.qubits) for op in ops)
        end
        @test AG.verify_encoder(code, encoder; basis=:zero)
        @test AG.verify_encoder(code, encoder; basis=:plus)
        for r in 0:T
            prefix = AG.AxisCode(r, Ly, Lz)
            prefix_encoder = AG.AxisEncoder(prefix, encoder.operations[1:7N*r], encoder.checkpoints[1:r])
            @test AG.verify_encoder(prefix, prefix_encoder; basis=:zero)
            @test AG.verify_encoder(prefix, prefix_encoder; basis=:plus)
        end
    end
end

@testset "signed verifier rejects corrupt circuits and schedules" begin
    code = AG.AxisCode(2, 3, 3)
    encoder = AG.axis_encoder(code)
    @test_throws ArgumentError AG.verify_encoder(code, encoder; basis=:unknown)
    @test !AG.verify_encoder(AG.AxisCode(1,3,3), encoder)
    @test !AG.verify_encoder(code, AG.AxisEncoder(code, encoder.operations, [62,126]))
    for index in (1, 19, 20, 63, 64, 82, 83, 126)
        ops = deepcopy(encoder.operations)
        op = ops[index]
        if op.gate === :H
            ops[index] = AG.TN.EncoderOperation(:H, [last(op.qubits)+1])
        else
            ops[index] = AG.TN.EncoderOperation(:CNOT, reverse(op.qubits))
        end
        mutated = AG.AxisEncoder(code, ops, encoder.checkpoints)
        @test !AG.verify_encoder(code, mutated; basis=:zero)
        @test !AG.verify_encoder(code, mutated; basis=:plus)
    end
    for replacement in (AG.TN.EncoderOperation(:Z, [1]),
                        AG.TN.EncoderOperation(:CNOT, [1,1]),
                        AG.TN.EncoderOperation(:CNOT, [1,0]),
                        AG.TN.EncoderOperation(:CNOT, [1,AG.nqubits(code)+1]))
        ops = copy(encoder.operations)
        ops[1] = replacement
        @test !AG.verify_encoder(code, AG.AxisEncoder(code,ops,encoder.checkpoints))
    end
    # A correctly sized fresh gate inserted into transfer 1 cannot touch transfer 2.
    ops = copy(encoder.operations)
    ops[19] = AG.TN.EncoderOperation(:H, [AG.edge(code,:x,2,0,0)])
    @test !AG.verify_encoder(code, AG.AxisEncoder(code,ops,encoder.checkpoints))
end
