using HydrodynamicStability, Test

linear = (du, u, p, t) -> (@. du = p * u)
linear_analytic = (u0, p, t) -> u0 * exp(p * t)

u0 = [1 / 2]
tspan = (0.0, 1.0)
p = 1.01
prob = ODEProblem(linear, u0, tspan, p)

sol = solve(prob, RK4(), dy=0.01)

analytic_u = linear_analytic.(u0, p, sol.y)
numerical_u = copy(analytic_u)
for i in eachindex(sol.u)
    numerical_u[i] = sol.u[i][1]
end

@test numerical_u ≈ analytic_u

