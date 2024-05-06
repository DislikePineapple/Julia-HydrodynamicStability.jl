@doc """
`Linear stability problem`

Define a linear stavility problem that is solved for a compressible flate boundary layer which statisfy the Blasis solution:

``Lϕ(y) = 0``

#Constructors


"""

@doc """
`Ordinary Differential Equations Problem`

Define the ODEs problem
"""
struct ODEProblem{F, uType, yType, P, K} <: AbstractODEProblem
    f::F
    u0::uType
    yspan::yType
    p::P
    kwarg::K

    function ODEProblem(f, u0, yspan, p = NullParameter(); kwargs...)
        new{typeof(f), typeof(u0), typeof(yspan), typeof(p), typeof(kwargs)}(
            f,
            u0,
            yspan,
            p,
            kwargs
        )
    end
end

"""
    `Boundary Value Problem`

Define a boundary Value problem

## Problem type

```julia
    BVProblem{F, BC, U0, Y, P, K} <: AbstractBVProblem
```

## Fields

- `f`: Function for the ordianry differential equations ``du = f(u,t,p)``.
- `bc`: Boundary conditions for the ODEs. Given as a function.
- `yspan`: Wall-normal direction span for the problem.
- `u0`: Initial conditions for the ODEs.
- `p`: The parameters for the problem. Defaults to `NullParameters`
- `kwargs`: The keyword arguments.

"""
struct BVProblem{F, BC, U0, Y, P, K} <: AbstractODEProblem
    f::F
    bc::BC
    u0::U0
    yspan::Y
    p::P
    kwargs::K

    function BVProblem(f, bc, u0, yspan, p = NullParameter(); kwargs...)
        new{typeof(f), typeof(bc), typeof(u0), typeof(yspan), typeof(p), typeof(kwargs)}(
            f,
            bc,
            u0,
            yspan,
            p,
            kwargs
        )
    end
end
