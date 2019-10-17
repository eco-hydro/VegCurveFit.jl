using Documenter, nlminb

makedocs(;
    modules=[nlminb],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://kongdd/kongdd/nlminb.jl/blob/{commit}{path}#L{line}",
    sitename="nlminb.jl",
    authors="Dongdong Kong",
    assets=String[],
)

deploydocs(;
    repo="kongdd/kongdd/nlminb.jl",
)
