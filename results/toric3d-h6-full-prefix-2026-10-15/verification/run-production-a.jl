include(joinpath(@__DIR__, "..", "run_study.jl"))
for arguments in (
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-x-L11-a", "--sizes", "11", "--rates", "0.035,0.04,0.045", "--shots", "1000", "--seed", "2026101573", "--batch-size", "50"],
)
    status = H6FullPrefixStudy.main(arguments)
    status == 0 || exit(status)
end
