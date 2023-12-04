using HydrodynamicStability, Test
import UnPack: @unpack

include("type.jl")
include("equations.jl")

# include("baseflow.jl")

tspan = (0, 15)
u₀ = [0, 0, 0.3] # the initial guess

bvp = BVProblem(blasius!, blasius_bc!, u₀, tspan)
sol = solve(bvp, Shooting(), dy = 0.01)

# @show(sol.u[1][3])
@test isapprox(sol.u[1][3], λ₀, atol = 1e-4)

Ma = 3
Te = 250
Re = Inf
fs = FreeStream(Ma, Te, Re)
u₀ = [0, 0, 0.3779, 0, 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2, 0]
tspan = (0.0, 15.0)

bvp1 = BVProblem(similarity!, similarity_bc!, u₀, tspan, fs)
sol1 = solve(bvp1, Shooting(), dy = 0.01)

# T₀ = 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2
# λᵤ = λ₀ * C(T₀, Te)^(-1 / 2)
# @show(λᵤ)
# @show(sol1.u[1][3])
@test isapprox(sol1.u[1][3], 0.3913, atol = 1e-4)
