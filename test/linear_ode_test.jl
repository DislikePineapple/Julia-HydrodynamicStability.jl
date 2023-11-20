using HydrodynamicStability, Test

linear = (du, u, p, t) -> (@. du = p * u)
linear_analytic = (u0, p, t) -> u0 * exp(p * t)

u0 = [1 / 2]
tspan = (0.0, 1.0)
p = 1.01
prob = ODEProblem(linear, u0, tspan, p)

sol = solve(prob, RK4(), dy = 0.01)

analytic_u = linear_analytic.(u0, p, sol.y)

@test [x[1] for x in sol.u] ≈ analytic_u
