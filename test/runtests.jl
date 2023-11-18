using Pkg
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "Core" || GROUP == "All"
    @time @safetestset "Aqua" include("aqua.jl")
    @time @safetestset "Linear ODE" include("linear_ode_test.jl")
end
