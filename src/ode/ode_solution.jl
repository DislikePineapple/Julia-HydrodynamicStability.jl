@doc """
Definde the ordinary differential equations solution
"""
struct ODESolution{uType, yType, P, A, uType2, RT} <: AbstractODESolution
    u::uType
    y::yType

    prob::P
    alg::A

    u_analytic::uType2
    retcode::RT
end

ODESolution(u, y) = ODESolution(u, y, nothing, nothing)
ODESolution(u, y, prob, alg) = ODESolution(u, y, prob, alg, nothing, nothing)

## Solution interface
@recipe function f(sol::AbstractODESolution)
    sol.y, u2matrix(sol)
end

function u2matrix(sol)
    u = zeros(eltype(eltype(sol.u)), length(sol.u), length(sol.u[1]))
    for i in eachindex(sol.u)
        u[i, :] = sol.u[i]
    end
    return u
end
