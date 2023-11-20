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
struct ODEProblem{F,uType,yType,P,K} <: AbstractODEProblem
    f::F
    u0::uType
    yspan::yType
    p::P
    kwarg::K

    function ODEProblem(f, u0, yspan, p = NullParameter(); kwargs...)
        new{typeof(f),typeof(u0),typeof(yspan),typeof(p),typeof(kwargs)}(
            f,
            u0,
            yspan,
            p,
            kwargs,
        )
    end
end

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
struct BVProblem{F,BC,Y,U0,P,K} <: AbstractBVProblem
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

    function BVProblem(f, bc, yspan, u0, p = NullParameter(); kwargs...)
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
