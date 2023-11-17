using Documenter, HydrodynamicStability

include("pages.jl")

makedocs(
    sitename = "HydrodynamicStability.jl",
    authors = "Sheng Yang",
    # modules = [HydrodynamicStability],
    # clean = true, doctest = false, linkcheck = true,
    # format=Documenter.HTML(
    # # assets=["assets/favicon.ico"],
    # canonical = "https://super-popar-bear.github.io/HydrodynamicStability/"),
    pages = pages,
)

deploydocs(
    repo = "github.com/Super-Popar-Bear/HydrodynamicStability.git";
    push_preview = true,
)
