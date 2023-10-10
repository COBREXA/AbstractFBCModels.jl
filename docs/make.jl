using Documenter, Literate
import AbstractFBCModels

function replace_listings(str)
    makelisting(line) = begin
        filename = line[8:end]
        open(joinpath("..", filename), "r") do f
            "```julia\n" * join(readlines(f) .* "\n") * "```\n"
        end
    end
    replace(str, r"^##LIST [^\s]*$"m => makelisting)
end

examples =
    sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/COBREXA/AbstractFBCModels.jl/blob/master",
        postprocess = replace_listings,
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

makedocs(
    modules = [AbstractFBCModels],
    clean = false,
    format = Documenter.HTML(
        ansicolor = true,
        canonical = "https://cobrexa.github.io/AbstractFBCModels.jl/stable/",
    ),
    sitename = "AbstractFBCModels.jl",
    linkcheck = false,
    pages = ["README" => "index.md"; example_mds; "Reference" => "reference.md"],
)

deploydocs(
    repo = "github.com/COBREXA/AbstractFBCModels.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
