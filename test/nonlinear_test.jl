using HydrodynamicStability, Test
f(t, p) = t .* t - p
df(u, t, p) = 2t
p = 2

prob1 = NonlinearProblem((f, df), 1; p = p)
sol1 = solve(prob1, Newton(); maxiters = 100, abstol = 1e-9)
@test sol1.t[1] ≈ √p

prob2 = NonlinearProblem(f, 1; p = p)
for alg in (Secant, Muller)
    sol2 = solve(prob2, alg(); maxiters = 100, abstol = 1e-9)
    @test sol2.t[1] ≈ √p
end

prob3 = NonlinearProblem(f, [1, 2]; p = p)
for alg in (Bisection, Falsi)
    sol3 = solve(prob3, alg(); abstol = 1e-9)
    @test sol3.t[1] ≈ √p
end
