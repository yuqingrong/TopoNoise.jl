include(joinpath(@__DIR__, "..", "run_study.jl"))
for arguments in (
    ["--channel", "x_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/pilot-x", "--sizes", "3,7,11", "--rates", "0,.01,.02,.035,.05,.065,.08", "--shots", "10", "--seed", "2026101561", "--batch-size", "10"],
    ["--channel", "z_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/pilot-z", "--sizes", "3,7,11", "--rates", "0,.01,.02,.03,.04,.05,.06", "--shots", "10", "--seed", "2026101562", "--batch-size", "10"],
    ["--channel", "z_only", "--output", "results/toric3d-h6-full-prefix-2026-10-15/pilot-z-refine", "--sizes", "3,7,11", "--rates", ".025,.03,.035,.04,.045,.05", "--shots", "100", "--seed", "2026101563", "--batch-size", "25"],
)
    status = H6FullPrefixStudy.main(arguments)
    status == 0 || exit(status)
end
