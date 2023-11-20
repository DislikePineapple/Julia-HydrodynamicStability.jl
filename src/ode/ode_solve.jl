# struct ODEcache{yType,P,Alg}
#     y::yType

#     prob::P
#     alg::Alg
# end

# solve(prob::ODEProblem, alg, args...; dy, kwargs...) = solve(init(prob, alg, args...; dy, kwargs...))
# solve(cache) = solve(cache, cache.alg)

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

# function solve(cache::ODEcache, alg::RK4)
#     @unpack f, u0, p = cache.prob
#     u = rk4(f, u0, cache.y, p)
#     sol = ODESolution(u, cache.y, prob, alg)
#     return sol
# end

# function solve(prob::BVProblem, alg::Shooting, args...; dy, kwargs...)
#     @unpack f, bc, yspan, u0, p = prob
#     @unpack ivp, iter = alg
#     @unpack tol, δ = iter

#     y = yspan[1]:dy:yspan[2]

#     # Blasius = shooting(f, bc, u0, y, p)
# end

# # shooting(fun, bcFun, u0, y) = shooting(fun, bcFun, u0, y, RK4())
# function shooting(prob:BVProblem, ::IvpAlgorithm; eps=1e-7, type=Float64, δ=0.0001)
#     function ivp(guess)
#         bc = zeros(type, length(guess))
#         f = method(fun, [ic; guess], t, type=type, para=para, para_fun=para_fun)
#         bcFun(bc, f, para_bc, t)
#         return bc
#     end
#     method(fun, [ic; secant(ivp, guess, eps=eps, δ=δ, type=type)], t, type=type, para=para, para_fun=para_fun)
#     # return IVP
# end
