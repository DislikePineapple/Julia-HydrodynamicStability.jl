"""
Define a partial differential equation problem used for the calculation of baseflow.
The first step is only consider the parabolishell equation, which is mainly used for the boundary layer problem.
```math
V_{xx}\frac{∂^2 u_{i,j}}{∂ x^2} + A\frac{∂ u_{i,j}}{∂ x} + (D + \frac{3\varGamma}{2Δt})u_{i,j} = \frac{\varGamma}{2Δt}(4u_{i-1,j}-u_{i-2,j}). 
```
where `u` is the baseflow, `x` is the spatial coordinate, `t` is the time, the parameter should take the form of [V_{xx},A,D,\varGamma], in where each element is a two-dimetional array.
"""

struct HeatProblem{G,I,FP,BCP,K} <: AbstractPDEProblem
    grid::G

    ic::I
    fun_para::FP
    bc_para::BCP

    kwarg::K

    function HeatProblem(grid, ic, fun_para, bc_para, kwarg...)
        new{
            typeof(grid),
            typeof(ic),
            typeof(fun_para),
            typeof(bc_para),
            typeof(kwarg),
        }(
            grid,
            ic,
            fun_para,
            bc_para,
            kwarg,
        )
    end
end
