using HydrodynamicStability, Test

function problem1!(du, u, t, p)
    du[1] = u[2]
    du[2] = u[1] - u[1]^2
end

function bc1!(residual, u, t, p)
    residual[1] = u[begin][1] - 1
    residual[2] = u[end][1] - 4
end

trange = (0, 1)
D, tspan = chebyshevshift(64, trange)

u₀ = [1.0, 0.0]
bvp = BVProblem(problem1!, bc1!, u₀, tspan)
ret = solve(bvp, Shooting())
sol = flatten_vector(ret.u)

fun1 = D * sol[1, :] - sol[2, :]
fun2 = D * sol[2, :] - (sol[1, :] - sol[1, :] .^ 2)
@test isapprox(maximum(abs, fun1), 0, atol = 1e-5)
@test isapprox(maximum(abs, fun2), 0, atol = 1e-5)

##* Blasius Equation ------------------------------------------------

function blasius!(du, u, t, p)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
end

function blasius_bc!(residual, u, t, p)
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[end][2] - 1
end

trange = (0, 15)
D, tspan = chebyshevshift(64, trange)

u₀ = [0, 0, 0.3] # the initial guess
bvp = BVProblem(blasius!, blasius_bc!, u₀, tspan)
ret = solve(bvp, Shooting())
sol = flatten_vector(ret.u)

fun1 = D * sol[1, :] - sol[2, :]
fun2 = D * sol[2, :] - sol[3, :]
fun3 = D * sol[3, :] + 1 / 2 * sol[1, :] .* sol[3, :]

@test isapprox(maximum(abs, fun1), 0, atol = 1e-4)
@test isapprox(maximum(abs, fun2), 0, atol = 1e-4)
@test isapprox(maximum(abs, fun3), 0, atol = 1e-4)

@test isapprox(sol[3, 1], 0.33206, atol = 1e-5)
