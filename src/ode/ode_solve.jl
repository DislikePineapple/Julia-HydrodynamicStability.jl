function solve(prob::ODEProblem, alg::RK4, args...; dy, kwargs...)
    @unpack f, u0, yspan, p = prob
    y = collect(yspan[1]:dy:yspan[2])

    u = Array{typeof(u0)}(undef, length(y))

    u[1] = u0

    s1 = similar(u0)
    s2 = similar(u0)
    s3 = similar(u0)
    s4 = similar(u0)

    for i = 1:length(y)-1
        h = y[i+1] - y[i]
        f(s1, u[i], p, y[i])
        f(s2, u[i] .+ h / 2 * s1, p, y[i] + h / 2)
        f(s3, u[i] .+ h / 2 * s2, p, y[i] + h / 2)
        f(s4, u[i] .+ h * s3, p, y[i] + h)
        u[i+1] = u[i] .+ h / 6 * (s1 .+ 2 * s2 .+ 2 * s3 .+ s4)
    end

    ODESolution(u, y)
end

function solve(
    prob::BVProblem,
    alg::Shooting,
    args...;
    dy,
    abstol = nothing,
    maxiters = 1000,
    kwargs...,
)
    @unpack f, bc, yspan, u0, p = prob
    @unpack ivp, iter = alg

    y = collect(yspan[1]:dy:yspan[2])
    bc!(residual, u) = bc(residual, u, p, y)

    function f_non(u0, p)
        residual = similar(u0)
        sol = solve(ODEProblem(f, u0, yspan, p), ivp, dy = dy)
        bc!(residual, sol.u)
        return residual
    end

    ic = solve(NonlinearProblem(f_non, u0, p), abstol = abstol, maxiters = maxiters)
    solve(ODEProblem(f, ic.t, yspan, p), ivp, dy = dy)
end
