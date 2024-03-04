"""
    solve(prob::PDEProblem, alg::FDM, args...; kwarg...)

Solve the PDE problem using the backward finite difference method with fixed step-length grid.
```math
\frac{∂^2u}{∂x^2} - λ\frac{∂u}{∂t} = 0 
```
"""
abstract type PDEAlgorithm <: AbstractAlgorithm end
struct CrankNicolson <: PDEAlgorithm end

function solve(prob::HeatProblem, alg::FDM, args...; kwarg...)
    @unpack grid, ic, fun_para, bc_para = prob
    _, tn, xn = size(grid)
    uType = eltype(ic)
    flow = zeros(uType, 2, tn, xn)
    p = Array{Array{uType}}(undef, xn)

    Δt = grid[1, 2, 1] - grid[1, 1, 1]

    for i = 1:tn
        if i == 1
            flow[1, i, :] = ic
            continue
        elseif i == 2
            for j = 1:xn
                p[j] = [
                    -fun_para[3][i, j] / fun_para[1][i, j] -
                    fun_para[4][i, j] / fun_para[1][i, j] / Δt,
                    fun_para[2][i, j] / fun_para[1][i, j],
                    flow[1, i-1, j] * fun_para[4][i, j] / fun_para[1][i, j] / Δt,
                ]
            end
        else
            for j = 1:xn
                p[j] = [
                    -fun_para[3][i, j] / fun_para[1][i, j] -
                    3 / 2 * fun_para[4][i, j] / fun_para[1][i, j] / Δt,
                    fun_para[2][i, j] / fun_para[1][i, j],
                    (4 * flow[1, i-1, j] - flow[1, i-2, j]) / 2 * fun_para[4][i, j] /
                    fun_para[1][i, j] / Δt,
                ]
            end
        end

        function heat_function!(M, F, p, t)
            M[1, 2] = 1
            M[2, 1] = p[1]
            M[2, 2] = p[2]

            F[1] = 0
            F[2] = p[3]
        end

        function heat_bc!(M0, Mend, u)
            M0[2, :] .= 0
            M0[2, 1] = 1
            Mend[2, :] .= 0
            Mend[2, end-1] = 1

            u[2] = bc_para[1][i]
            u[end] = bc_para[2][i]
        end

        ode_prob = BVProblem(heat_function!, heat_bc!, [0.0, 0.0], grid[2, i, :], p)
        sol = solve(ode_prob, FDM())
        for j = 1:xn
            flow[1, i, j] = sol.u[j][1]
            flow[2, i, j] = sol.u[j][2]
        end
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function solve(prob::HeatProblem, alg::CrankNicolson, args...; kwarg...)
    @unpack grid, ic, fun_para, bc_para = prob
    _, tn, xn = size(grid)
    uType = eltype(ic)
    flow = zeros(uType, 2, tn, xn)
    p = Array{Array{uType}}(undef, xn)

    Δt = grid[1, 2, 1] - grid[1, 1, 1]

    for i = 1:tn
        if i == 1
            flow[1, i, :] = ic
            continue
        end
        udoubleprime = central_difference(
            central_difference(flow[1, i-1, :], grid[2, i-1, :]),
            grid[2, i-1, :],
        )
        if i == 2
            for j = 1:xn
                p[j] = [
                    -2 * fun_para[3][i, j] / fun_para[1][i, j] -
                    2 * fun_para[4][i, j] / fun_para[1][i, j] / Δt,
                    2 * fun_para[2][i, j] / fun_para[1][i, j],
                    2 * flow[1, i-1, j] * fun_para[4][i, j] / fun_para[1][i, j] / Δt -
                    udoubleprime[j],
                ]
            end
        else
            for j = 1:xn
                p[j] = [
                    -2 * fun_para[3][i, j] / fun_para[1][i, j] -
                    3 * fun_para[4][i, j] / fun_para[1][i, j] / Δt,
                    2 * fun_para[2][i, j] / fun_para[1][i, j],
                    (4 * flow[1, i-1, j] - flow[1, i-2, j]) * fun_para[4][i, j] /
                    fun_para[1][i, j] / Δt - udoubleprime[j],
                ]
            end
        end

        function heat_function!(M, F, p, t)
            M[1, 2] = 1
            M[2, 1] = p[1]
            M[2, 2] = p[2]

            F[1] = 0
            F[2] = p[3]
        end

        function heat_bc!(M0, Mend, u)
            M0[2, :] .= 0
            M0[2, 1] = 1
            Mend[2, :] .= 0
            Mend[2, end-1] = 1

            u[2] = bc_para[1][i]
            u[end] = bc_para[2][i]
        end

        ode_prob = BVProblem(heat_function!, heat_bc!, [0.0, 0.0], grid[2, i, :], p)
        sol = solve(ode_prob, FDM())
        for j = 1:xn
            flow[1, i, j] = sol.u[j][1]
            flow[2, i, j] = sol.u[j][2]
        end
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end