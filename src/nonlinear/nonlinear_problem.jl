"""
    `NonlinearProblem`

Defines the nonlinear problem.

## Mathematical Specification of the Nonlinear Problem

Problem of nonlinear equation, if
```math
    f(x_0) = 0,
```
function ``f(x)`` has a root when ``x=x_0``.

## Problem type

```julia
    NonlinearProblem{F, DF, tType, P, K} <: AbstractProblem
```

## Fields

- `f`: Function for the nonlinear equations.
- `t0`: Initial guess for the nonlinear equations.
- `p`: The parameters for the problem. Defaults to `NullParameters`
- `kwargs`: The keyword arguments.

"""
struct NonlinearProblem{F, tType, P, K} <: AbstractProblem
    f::F
    t0::tType

    p::P
    kwarg::K

    function NonlinearProblem(f, t0; p = NullParameter(), kwargs...)
        new{typeof(f), typeof(t0), typeof(p), typeof(kwargs)}(f, t0, p, kwargs)
    end
end
