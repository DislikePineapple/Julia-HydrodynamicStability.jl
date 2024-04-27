"""
Define the solution of the partial differential equation.
"""
struct HeatSolution{uType, yType, P, A, uType2, RT} <: AbstractPDESolution
    grid::yType
    flow::uType

    prob::P
    alg::A

    u_analytic::uType2
    retcode::RT
end

HeatSolution(grid, flow) = HeatSolution(grid, flow, nothing, nothing)
HeatSolution(grid, flow, prob, alg) = HeatSolution(grid, flow, prob, alg, nothing, nothing)
