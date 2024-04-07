"""
    solve(prob::PDEProblem, alg::FDM, args...; kwarg...)

Solve the PDE problem using the backward finite difference method with fixed step-length grid.
```math
\frac{∂^2u}{∂x^2} - λ\frac{∂u}{∂t} = 0 
```
"""
abstract type PDEAlgorithm <: AbstractAlgorithm end
# struct FDM <: PDEAlgorithm end

function solve(prob::HeatProblem, alg::FDM, args...; kwarg...)
    @unpack grid, ic = prob
    _, nt, nx = size(grid)
    Type = eltype(eltype(ic))
    O = length(ic[1])


    flow = zeros(Type, O, nt, nx)

    for j = 1:nt
        if j == 1
            for k = 1:nx
                flow[:, j, k] = ic[k]
            end
            continue
        else
            M, B = jac(prob, flow, j)
            B = M \ B
            for k = 1:nx
                flow[:, j, k] = B[(k-1)*O+1:k*O]
            end
        end
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function jac(prob, flow, tn)
    @unpack fun, fun_para, bc, bc_para, grid, ic = prob
    nx = size(grid)[3]
    T = eltype(eltype(ic))
    O = length(ic[1])

    M = zeros(T, nx * O, nx * O)
    B = zeros(T, nx * O)

    Vxx, A, Γ, D = [zeros(T, O, O) for _ = 1:4]
    F = zeros(T, O)

    x = grid[2, tn, :]
    ΔT = grid[1, tn, 1] - grid[1, tn-1, 1]

    tn == 2 ? C0 = 1 : C0 = 3 / 2

    for i in eachindex(x)
        position, dx = fun_position(x, i)
        for j in eachindex(x)
            if fun_para isa NullParameter || fun_para[1] isa Number
                fun(Vxx, A, Γ, D, F, fun_para)
            else
                fun(Vxx, A, Γ, D, F, fun_para[tn, j])
            end

            if i == j
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] +=
                    C2[position, 3] .* Vxx ./ dx^2 +
                    C1[position, 3] .* A ./ dx +
                    D +
                    C0 .* Γ ./ ΔT
                if tn == 2
                    B[(i-1)*O+1:i*O] = Γ * flow[:, tn-1, i] ./ ΔT + F
                else
                    B[(i-1)*O+1:i*O] =
                        Γ * (4 .* flow[:, tn-1, i] - flow[:, tn-2, i]) ./ 2 ./ ΔT + F
                end
            elseif j == i - 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] +=
                    C2[position, 1] .* Vxx ./ dx^2 + C1[position, 1] .* A ./ dx
            elseif j == i - 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] +=
                    C2[position, 2] .* Vxx ./ dx^2 + C1[position, 2] .* A ./ dx
            elseif j == i + 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] +=
                    C2[position, 4] .* Vxx ./ dx^2 + C1[position, 4] .* A ./ dx
            elseif j == i + 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] +=
                    C2[position, 5] .* Vxx ./ dx^2 + C1[position, 5] .* A ./ dx
            end
        end
    end

    M0 = M[1:O, :]
    Mend = M[end-O+1:end, :]

    bc(M0, Mend, B, bc_para[tn])

    M[1:O, :] = M0
    M[end-O+1:end, :] = Mend

    return M, B
end
