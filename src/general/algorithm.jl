abstract type NonlinearAlgorithm <: AbstractAlgorithm end
struct Bisection <: NonlinearAlgorithm end
struct Muller{T,Δ} <: NonlinearAlgorithm
    tol::T
    δ::Δ
end
Muller() = Muller(1e-9, 1e-7 + 1e-5im)
struct Secant{T,Δ} <: NonlinearAlgorithm
    tol::T
    δ::Δ
end

abstract type ODEsAlgorithm <: AbstractAlgorithm end
struct RK4 <: ODEsAlgorithm end

function rk4(f::Function, u0::AbstractArray, t::AbstractArray, p)
    u = zeros(typeof(u0), length(u0), length(t))
    s1 = zeros(typeof(u0), length(u0))
    s2 = zeros(typeof(u0), length(u0))
    s3 = zeros(typeof(u0), length(u0))
    s4 = zeros(typeof(u0), length(u0))

    u[:, 1] = u0 # initial conditions

    for i = 1:length(t)-1
        h = t[i+1] - t[i]
        f(s1, u[:, i], p, t[i])
        f(s2, u[:, i] .+ h / 2 * s1, p, t[i] + h / 2)
        f(s3, u[:, i] .+ h / 2 * s2, p, t[i] + h / 2)
        f(s4, u[:, i] .+ h * s3, p, t[i] + h)
        u[:, i+1] = u[:, i] .+ h / 6 * (s1 .+ 2 * s2 .+ 2 * s3 .+ s4)
    end
    return u
end

abstract type BVProblemAlgorithm <: AbstractAlgorithm end
struct Shooting{I,IT} <: BVProblemAlgorithm
    ivp::I
    iter::IT
end

Shooting() = Shooting(RK4(), Muller())
