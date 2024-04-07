"""
Define a partial differential equation problem used for the calculation of baseflow.
The first step is only consider the parabolishell equation, which is mainly used for the boundary layer problem.
```math
V_{xx}\frac{∂^2 u_{i,j}}{∂ x^2} + A\frac{∂ u_{i,j}}{∂ x} + D u_{i,j} + \frac{3Γ}{2Δt}u_{i,j} = \frac{Γ}{2Δt}(4u_{i-1,j}-u_{i-2,j}) + F. 
```
where `u` is the baseflow, `x` is the spatial coordinate, `t` is the time, the parameter should take the form of [V_{xx}, A, D, Γ, F], in where each element is a two-dimetional array.
"""

struct HeatProblem{G,F,BC,I,FP,BCP,K} <: AbstractPDEProblem
    grid::G

    fun::F
    bc::BC

    ic::I
    fun_para::FP
    bc_para::BCP

    kwarg::K

    function HeatProblem(grid, fun, bc, ic, fun_para, bc_para, kwarg...)
        new{
            typeof(grid),
            typeof(fun),
            typeof(bc),
            typeof(ic),
            typeof(fun_para),
            typeof(bc_para),
            typeof(kwarg),
        }(
            grid,
            fun,
            bc,
            ic,
            fun_para,
            bc_para,
            kwarg,
        )
    end
end
