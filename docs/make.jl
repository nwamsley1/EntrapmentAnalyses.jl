using EntrapmentAnalyses
using Documenter

DocMeta.setdocmeta!(EntrapmentAnalyses, :DocTestSetup, :(using EntrapmentAnalyses); recursive=true)

makedocs(;
    modules=[EntrapmentAnalyses],
    authors="Nathan Wamsley and contributors",
    sitename="EntrapmentAnalyses.jl",
    format=Documenter.HTML(;
        canonical="https://nathanwamsley.github.io/EntrapmentAnalyses.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => [
            "Getting Started" => "guide.md",
            "Tutorial" => "tutorial.md",
        ],
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(;
    repo="github.com/nathanwamsley/EntrapmentAnalyses.jl",
    devbranch="master",
)