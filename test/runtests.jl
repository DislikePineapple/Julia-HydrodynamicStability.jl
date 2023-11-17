using Pkg
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "Core" || GROUP == "All"
    @time @safetestset "Aqua" include("aqua.jl")
end
