@doc """
`Linear stability problem`

Define a linear stavility problem that is solved for a compressible flate boundary layer which statisfy the Blasis solution:

``Lϕ(y) = 0``

#Constructors


"""

@doc """
`Boundary Value Problem`

Define a boundary Value problem

### Fields

- `f`: Function for the ordianry differential equations Du=f(u,t,p).
- `bc`: Boundary conditions for the ODEs. Given as a function.
- `yspan`: Wall-normal direction span for the problem.
- `u0`: Initial conditions for the ODEs.
- `p`: The parameters for the problem. Defaults to `NullParameters`
- `kwargs`: The keyword arguments.
"""
struct BVProblem{isinplace,F,BC,Y,U0,P,K}
    <:AbstractBVProblem
    "f: Function for the ordianry differential equations Du=f(u,t,p)."
    f::F
    "bc: Boundary conditions for the ODEs. Given as a function."
    bc::BC
    "yspan: Wall-normal direction span for the problem."
    yspan::Y
    "u0: Initial conditions for the ODEs."
    u0::U0
    "p: Parameter for problem."
    p::P
    "kwargs: The keyword arguments."
    kwargs::K

    function BVProblem(f, bc, yspan, u0, p = NullParameters(); kwargs...)
        new{typeof(f),typeof(bc),typeof(yspan),typeof(u0),typeof(p),typeof(kwargs)}(
            f,
            bc,
            yspan,
            u0,
            p,
            kwargs,
        )
    end
end
