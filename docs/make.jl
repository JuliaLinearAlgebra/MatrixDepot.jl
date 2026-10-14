# use the index snapshot in test/data, as test/runtests.jl does, so the build needs no network
ENV["MATRIXDEPOT_URL_REDIRECT"] = "1"
ENV["MATRIXDEPOT_DATA"] = mktempdir()

using Documenter, MatrixDepot

DocMeta.setdocmeta!(MatrixDepot, :DocTestSetup, :(using MatrixDepot); recursive=true)

makedocs(
    sitename = "MatrixDepot.jl",
    modules = [MatrixDepot],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://JuliaLinearAlgebra.github.io/MatrixDepot.jl/stable/",
    ),
    checkdocs = :exports, # internal helpers carry docstrings not meant for the manual
    pages = [
        "Home" => "index.md",
        "Basics" => "basics.md",
        "Remote Collections" => "remote.md",
        "Matrices" => "matrices.md",
        "Regularization Problems" => "regu.md",
        "Random Graphs" => "graphs.md",
        "Groups" => "groups.md",
        "User Generators" => "user.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
)

deploydocs(repo = "github.com/JuliaLinearAlgebra/MatrixDepot.jl.git", push_preview = true)
