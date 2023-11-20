using HydrodynamicStability, Test
f(t, p) = t .* t - p
p = 2
prob1 = NonlinearProblem(f, 1, p)

for alg in (Secant, Muller,)
    sol = solve(prob1, alg(), tol=1e-9)
    @test sol.t[1] ≈ √p
end

prob2 = NonlinearProblem(f, [1, 2], p)
for alg in (Bisection, Falsi,)
    sol = solve(prob2, alg(), tol=1e-9)
    @test sol.t[1] ≈ √p
end
