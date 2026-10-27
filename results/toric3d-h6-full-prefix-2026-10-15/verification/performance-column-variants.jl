using Random, Statistics, Printf
function with_view(matrix::BitMatrix, frame::BitVector)
    result = falses(size(matrix, 1))
    for j in findall(frame)
        result .⊻= view(matrix, :, j)
    end
    result
end
function with_copy(matrix::BitMatrix, frame::BitVector)
    result = falses(size(matrix, 1))
    for j in findall(frame)
        result .⊻= matrix[:, j]
    end
    result
end
function timed(f, matrix, frame)
    f(matrix, frame)
    samples = [@elapsed f(matrix, frame) for _ in 1:9]
    return median(samples), @allocated(f(matrix, frame))
end
function main()
    rng = MersenneTwister(20261015)
    matrix = bitrand(rng, 3993, 4356)
    println("Random bit matrix 3993x4356; median of 9 hot calls")
    println("density,weight,equal,view_ms,copy_ms,speedup,view_bytes,copy_bytes")
    for p in (.05, .5)
        frame = BitVector(rand(rng, 4356) .< p)
        @assert with_view(matrix, frame) == with_copy(matrix, frame)
        view_time, view_bytes = timed(with_view, matrix, frame)
        copy_time, copy_bytes = timed(with_copy, matrix, frame)
        @printf("%.2f,%d,true,%.6f,%.6f,%.3f,%d,%d\n", p, count(frame),
            1000view_time, 1000copy_time, view_time/copy_time, view_bytes, copy_bytes)
    end
    println("COMPLETE")
end
main()
