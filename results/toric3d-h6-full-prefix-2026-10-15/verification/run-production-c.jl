include(joinpath(@__DIR__, "..", "run_study.jl"))
for arguments in (
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-x-small", "--sizes", "3,5,7", "--rates", "0.035,0.04,0.045,0.05,0.055,0.06,0.065", "--shots", "1000", "--seed", "2026101571", "--batch-size", "50"],
    ["--channel", "z_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-z", "--sizes", "3,5,7,9,11", "--rates", "0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04", "--shots", "1000", "--seed", "2026101576", "--batch-size", "50"],
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/production-x-L11-c", "--sizes", "11", "--rates", "0.06,0.065", "--shots", "1000", "--seed", "2026101575", "--batch-size", "50"],
)
    status = H6FullPrefixStudy.main(arguments)
    status == 0 || exit(status)
end
