include(joinpath(@__DIR__, "..", "run_study.jl"))
for arguments in (
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-x-L11-b", "--sizes", "11", "--rates", "0.05,0.055", "--shots", "1000", "--seed", "2026101574", "--batch-size", "50"],
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-x-L9", "--sizes", "9", "--rates", "0.035,0.04,0.045,0.05,0.055,0.06,0.065", "--shots", "1000", "--seed", "2026101572", "--batch-size", "50"],
)
    status = H6FullPrefixStudy.main(arguments)
    status == 0 || exit(status)
end
