using HydrodynamicStability, Test

function blasius!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
end

function blasius_bc!(residual, u, p, t)
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[end][2] - 1
end

tspan = (0, 15)
u₀ = [0, 0, 0.3] # the initial guess

bvp = BVProblem(blasius!, blasius_bc!, u₀, tspan)
sol = solve(bvp, Shooting(), dy = 0.01)

@test isapprox(sol.u[1][3], 0.3321, atol = 1e-4)
