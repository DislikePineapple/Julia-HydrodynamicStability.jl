"""
Define a partial differential equation problem used for the calculation of baseflow.
The first step is only consider the parabolishell equation, which is mainly used for the boundary layer problem.
```math
V_{xx}\frac{∂^2 u_{i,j}}{∂ x^2} + A\frac{∂ u_{i,j}}{∂ x} + D u_{i,j} + \frac{3Γ}{2Δt}u_{i,j} = \frac{Γ}{2Δt}(4u_{i-1,j}-u_{i-2,j}) + F. 
```
where `u` is the baseflow, `x` is the spatial coordinate, `t` is the time, the parameter should take the form of [V_{xx}, A, D, Γ, F], in where each element is a two-dimetional array.
"""

struct HeatProblem{G,F,I,IC,BC,FP,IP,BP,K} <: AbstractPDEProblem
    grid::G

    fun::F
    imhomo::I

    ic::IC
    bc::BC

    fun_para::FP
    im_para::IP
    bc_para::BP

    kwarg::K

    function HeatProblem(
        g,
        f,
        im,
        ic,
        bc,
        fp = NullParameter(),
        ip = NullParameter(),
        bp = NullParameter(),
        kwarg...,
    )
        new{
            typeof(g),
            typeof(f),
            typeof(im),
            typeof(ic),
            typeof(bc),
            typeof(fp),
            typeof(ip),
            typeof(bp),
            typeof(kwarg),
        }(
            g,
            f,
            im,
            ic,
            bc,
            fp,
            ip,
            bp,
            kwarg,
        )
    end
end

function HeatProblem(grid, funs, paras)
    fun, imhomo, bc = funs
    fun_para, im_para, ic, bc_para = paras
    HeatProblem(grid, fun, imhomo, ic, bc, fun_para, im_para, bc_para)
end

struct NSHeatProblem{G,Int,FUN,IM,IC,BC,FP,IMP,BCP,K} <: AbstractPDEProblem
    grid::G
    Ny::Int

    fun::FUN
    imhomo::IM

    ic::IC
    bc::BC

    fun_para::FP
    imhomo_para::IMP
    bc_para::BCP


    kwarg::K

    function NSHeatProblem(
        grid,
        Ny,
        fun,
        imhomo,
        ic,
        bc,
        fun_para,
        im_para,
        bc_para,
        kwarg...,
    )
        new{
            typeof(grid),
            Int,
            typeof(fun),
            typeof(imhomo),
            typeof(ic),
            typeof(bc),
            typeof(fun_para),
            typeof(im_para),
            typeof(bc_para),
            typeof(kwarg),
        }(
            grid,
            Ny,
            fun,
            imhomo,
            ic,
            bc,
            fun_para,
            im_para,
            bc_para,
            kwarg,
        )
    end
end

function NSHeatProblem(grid, Ny::Int, funs, paras)
    fun, imhomo, bc = funs
    fun_para, im_para, ic, bc_para = paras
    NSHeatProblem(grid, Ny, fun, imhomo, ic, bc, fun_para, im_para, bc_para)
end