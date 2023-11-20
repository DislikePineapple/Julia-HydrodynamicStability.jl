@doc """
Definde the ordinary differential equations solution

"""
struct ODESolution{uType,yType,P,A,uType2,RT} <: AbstractODESolution
    u::uType
    y::yType

    prob::P
    alg::A

    u_analytic::uType2
    retcode::RT
end

ODESolution(u, y) = ODESolution(u, y, nothing, nothing)
ODESolution(u, y, prob, alg) = ODESolution(u, y, prob, alg, nothing, nothing)

@recipe function f(sol::AbstractODESolution)
    sol.y, [x[1] for x in sol.u]
end
