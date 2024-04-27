using Documenter, HydrodynamicStability

include("pages.jl")

makedocs(
    sitename = "HydrodynamicStability.jl",
    authors = "Sheng Yang",
    pages = pages
)

deploydocs(
    repo = "github.com/DislikePineapple/Julia-HydrodynamicStability.jl.git";
    push_preview = true
)
