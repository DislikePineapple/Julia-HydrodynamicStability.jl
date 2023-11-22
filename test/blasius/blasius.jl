# include("BaseFlow.jl")
include("constant.jl")

using HydrodynamicStability, Test

function blasius!(du, u, p, t) # u[1] = f, u[2] = f', u[3] = f'' incompressible blasius
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
end

function blasius_bc!(residual, u, p, t)
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[end][2] - 1 # f' → 1 as y → ∞
end

tspan = (0, 15)
u₀ = [0, 0, 0.3] # the initial guess

bvp = BVProblem(blasius!, blasius_bc!, u₀, tspan)
sol = solve(bvp, Shooting(), dy = 0.1)

@test isapprox(sol.u[1][3], λ_0, atol = 1e-4)

C(T, Te) = T^(1 / 2) * (1 + Se / Te) / (T + Se / Te)
Cprime_over_C(T, Tprime, Te) = 1 / 2 * Tprime / T * (Se / Te - T) / (T + Se / Te)
D(T, Te) = C(T, Te) * (1 / 2 / T - 1 / (T + Se / Te))
Dprime(T, Tprime, Te) =
    -Tprime * (Se + Te * T) * (Se^2 + 6 * Se * Te * T - 3 * Te^2 * T^2) /
    (4 * T^(3 / 2) * (Se + Te * T)^3)

function similarity!(du, u, p, t) # u[1] = f, u[2] = f', u[3] = g, u[4] = f'', u[5] = g'; f' = U, g = T , p[1] = Ma, p[2] = Te
    du[1] = u[2]
    du[2] = u[4]
    du[3] = u[5]
    du[4] = -1 / 2 / C(u[3], p[2]) * u[1] * u[4] - Cprime_over_C(u[3], u[5], p[2]) * u[4]
    du[5] =
        -Pr / 2 / C(u[3], p[2]) * u[1] * u[5] - Cprime_over_C(u[3], u[5], p[2]) * u[5] -
        Pr * (Γ - 1) * p[1]^2 * u[4]^2
end

function similarity_bc!(residual, u, p, t)
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[begin][3] - (1 + sqrt(Pr) * (Γ - 1) * p[1]^2 / 2)
    residual[4] = u[end][2] - 1 # U → 1 as y → ∞
    residual[5] = u[end][3] - 1 # T → 1 as y → ∞
end

Ma = 3
Te = 250

u₀ = [0, 0, 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2, 0.33, 0]
tspan = (0.0, 15.0)
p = [Ma, Te]

bvp1 = BVProblem(similarity!, similarity_bc!, u₀, tspan, p)
sol1 = solve(bvp1, Shooting(), dy = 0.1)

@test isapprox(sol1.u[1][4], 0.3913, atol = 1e-4)
