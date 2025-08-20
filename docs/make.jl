using DiracBilinears
using Documenter

DocMeta.setdocmeta!(DiracBilinears, :DocTestSetup, :(using DiracBilinears); recursive=true)

makedocs(;
    modules=[DiracBilinears],
    authors="Tatsuya Miki <tatsuya.miki.1201@gmail.com> and contributors",
    sitename="DiracBilinears.jl",
    format=Documenter.HTML(;
        canonical="https://TatsuyaMiki.github.io/DiracBilinears.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "License and Citation" => "license_citation.md"
    ],
)

deploydocs(;
    repo="github.com/TatsuyaMiki/DiracBilinears.jl",
    devbranch="main",
)
