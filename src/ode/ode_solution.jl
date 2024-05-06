"""
    `ODE Solution`

Definde the ordinary differential equations solution

## Solution type

```julia
    ODESolution{uType, yType, P, A, uType2, RT} <: AbstractODESolution
```

## Fields

- `u`: The solution of the ODEs.
- `y`: The wall-normal direction span for the problem.
- `prob`: The problem that was solved.
- `alg`: The algorithm that was used to solve the problem.
- `u_analytic`: The analytic solution of the ODEs.
- `retcode`: The return code of the solver.

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
