using TopoNoise
using Documenter

DocMeta.setdocmeta!(TopoNoise, :DocTestSetup, :(using TopoNoise); recursive=true)

makedocs(;
    modules=[TopoNoise],
    authors="yuqingrong",
    sitename="TopoNoise.jl",
    format=Documenter.HTML(;
        canonical="https://yuqingrong.github.io/TopoNoise.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Noisy trajectories" => "trajectories.md",
    ],
)

deploydocs(;
    repo="github.com/yuqingrong/TopoNoise.jl",
    devbranch="main",
)
