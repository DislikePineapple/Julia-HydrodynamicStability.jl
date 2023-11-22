abstract type NonlinearAlgorithm <: AbstractAlgorithm end
struct Bisection <: NonlinearAlgorithm end
struct Falsi <: NonlinearAlgorithm end
struct Muller <: NonlinearAlgorithm end
struct Secant <: NonlinearAlgorithm end

solve(prob::NonlinearProblem; kwarg...) = solve(prob::NonlinearProblem, Secant(); kwarg...)

function solve(
    prob::NonlinearProblem,
    alg::Secant,
    arg...;
    abstol = nothing,
    reltol = nothing,
    maxiters = 1000,
    kwarg...,
)
    t = float(prob.t0)
    f(t) = prob.f(t, prob.p)
    T = typeof(t)

    atol =
        abstol !== nothing ? abstol :
        real(oneunit(eltype(T))) * (eps(real(one(eltype(T)))))^(4 // 5)
    rtol = reltol !== nothing ? reltol : eps(real(one(eltype(T))))^(4 // 5)

    if t isa Number
        to = oftype(one(eltype(t)), Inf)
    else
        to = map(t -> oftype(one(eltype(t)), Inf), t)
    end

    for _ = 1:maxiters
        if t isa Number
            u = f(t)
            du = ForwardDiff.derivative(f, t)
        elseif t isa AbstractArray
            u = f(t)
            du = ForwardDiff.jacobian(f, t)
            # du = Jacobin(f, u, t, 1e-5)
            # du = FiniteDiff.finite_difference_jacobian(f, t)
        else
            error("Secant only supports Number and AbstactVector types.")
        end

        iszero(u) && return NonlinearSolution(t, prob, alg)

        Δt = du \ u
        t -= Δt

        isapprox(t, to, atol = atol, rtol = rtol) && return NonlinearSolution(t, prob, alg)
        to = t
    end
    error("Failed to converge in $maxiters iteu_righttions, and t = $t")
end

# function Jacobin(f, u0, t0, δ)
#     N = length(t0)
#     J = zeros(eltype(t0), N, N)
#     @show(N)

#     for i in eachindex(t0)
#         t = copy(t0)
#         t[i] = t[i] + δ
#         u = f(t)
#         @show(i)
#         @show(t)
#         @show(u)
#         J[i, :] = (u - u0) ./ δ
#     end
#     # @show(J)
#     return J
# end

function solve(
    prob::NonlinearProblem,
    alg::Muller,
    arg...;
    abstol = 1e-9,
    maxiters = 1000,
    δ = 1e-5,
    kwargs...,
)
    @unpack f, t0, p = prob
    t0 = float(t0)

    for _ = 1:maxiters
        t = [t0 - δ t0 + δ t0]
        u = f.(t, p)

        f12 = (u[1] - u[2]) / (t[1] - t[2])
        f13 = (u[1] - u[3]) / (t[1] - t[3])
        f23 = (u[2] - u[3]) / (t[2] - t[3])

        f123 = (f23 - f12) / (t[1] - t[3])

        w = f12 + f13 - f23
        sqrt = w^2 - 4 * u[3] / f123

        abs(w - sqrt) < abs(w + sqrt) ? denoms = w + sqrt : denoms = w - sqrt

        t0 = t0 - 2 * u[3] / denoms
        sum(abs, u[3]) < abstol && return NonlinearSolution(t0, prob, alg)
    end
    error("Failed to converge in $maxiters iteu_righttions, and t = $t0")
end

function solve(prob::NonlinearProblem, alg::Bisection, arg...; abstol = 1e-9, kwarg...)
    @unpack f, t0, p = prob
    t0 = float(t0)
    left = t0[1]
    right = t0[2]

    u_left = f(left, p)
    u_right = f(right, p)
    u_left .* u_right > 0 && error("Space not include zeros")

    sol = copy(left)
    while (right - left) / 2 > abstol / 100
        sol = (left + right) / 2
        u_sol = f(sol, p)
        abs(u_sol) < abstol && return NonlinearSolution(sol, prob, alg)
        u_left * u_sol < 0 ? (right = sol; u_right = u_sol) : (left = sol; u_left = u_sol)
    end
    error("Not find zeros, and sol = $sol")
end

function solve(prob::NonlinearProblem, alg::Falsi, arg...; abstol = 1e-9, kwarg...)
    @unpack f, t0, p = prob
    t0 = float(t0)
    left = t0[1]
    right = t0[2]

    u_left = f(left, p)
    u_right = f(right, p)
    u_left .* u_right > 0 && error("Space not include zeros")

    sol = copy(left)
    while (right - left) > abstol / 100
        add = (right - left) .* u_left ./ (u_left - u_right)
        sol = left + add
        u_so = f(sol, p)
        abs(u_so) < abstol && return NonlinearSolution(sol, prob, alg)
        u_left * u_so < 0 ? (right = sol; u_right = u_so) : (left = sol; u_left = u_so)
    end
    error("Not find zeros, and sol = $sol")
end
