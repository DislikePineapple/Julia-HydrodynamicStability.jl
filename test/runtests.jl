using Pkg
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

@time begin
    if GROUP == "Core" || GROUP == "All"
        @time @safetestset "Aqua" include("aqua.jl")
        @time @safetestset "Nonlinear" include("nonlinear_test.jl")
        @time @safetestset "Linear ODE" include("linear_ode_test.jl")
        @time @safetestset "Blasius Solution" include("blasius/blasius.jl")
        @time @safetestset "Eigenvalue" include("heating/heating.jl")
    end
end
