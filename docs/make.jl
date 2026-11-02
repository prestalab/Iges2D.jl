using Iges2D
using Documenter

DocMeta.setdocmeta!(Iges2D, :DocTestSetup, :(using Iges2D); recursive = true)

makedocs(;
    modules = [Iges2D],
    authors = "PrestaLab <prestalabx@gmail.com> and contributors",
    sitename = "Iges2D.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://prestalab.github.io/Iges2D.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo = "github.com/prestalab/Iges2D.jl",
    devbranch = "main",
)
