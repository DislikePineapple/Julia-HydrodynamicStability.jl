using HydrodynamicStability, Test

linear = (du, u, t, p) -> (@. du = p * u)
linear_analytic = (u0, t, p) -> u0 * exp(p * t)

u0 = [1 / 2]
tspan = range(0.0, 1.0, 101)

p = 1.01
prob = ODEProblem(linear, u0, tspan, p)

sol = solve(prob, RK4())

analytic_u = linear_analytic.(u0, sol.y, p)

@test [x[1] for x in sol.u] ≈ analytic_u
