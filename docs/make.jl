using Documenter
import AbstractFBCModels

makedocs(
    modules = [AbstractFBCModels],
    clean = false,
    format = Documenter.HTML(
        ansicolor = true,
        canonical = "https://cobrexa.github.io/AbstractFBCModels.jl/stable/",
    ),
    sitename = "AbstractFBCModels.jl",
    linkcheck = false,
    pages = ["README" => "index.md"; "Reference" => "reference.md"],
    strict = [:missing_docs, :cross_references, :example_block],
)

deploydocs(
    repo = "github.com/COBREXA/AbstractFBCModels.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
