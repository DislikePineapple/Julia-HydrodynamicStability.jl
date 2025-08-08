"""
    solve(prob::PDEProblem, alg::FDM, args...; kwarg...)

Solve the PDE problem using the backward finite difference method with fixed step-length grid.
```math
\frac{∂^2u}{∂x^2} - λ\frac{∂u}{∂t} = 0 
```
"""

# using FiniteDiff to solve the parabolic PDE problem
abstract type PDEAlgorithm <: AbstractPDEAlgorithm end
abstract type NPDEAlgorithm <: PDEAlgorithm end

struct PFDM <: PDEAlgorithm end
struct PSFDM <: PDEAlgorithm end
struct NPSFDM <: PDEAlgorithm end

function solve(prob::HeatProblem, alg::PDEAlgorithm, args...; kwarg...)
    @unpack grid, ic = prob
    _, nt, nx = size(grid)
    Type = eltype(ic)
    nv = size(ic)[1]

    flow = zeros(Type, nv, nt, nx)

    for j in 1:nt
        if j == 1
            flow[:, j, :] = ic
            continue
        else
            M, B = jac(prob, alg, flow, j)
            B = M \ B
            for k in 1:nx
                rk = (1:nv) .+ ((k - 1) * nv)
                flow[:, j, k] = B[rk]
            end
        end
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function solve(
        prob::NSHeatProblem,
        alg::NPSFDM,
        args...;
        abstol = nothing,
        reltol = nothing,
        iterate = false,
        showprogress = false,
        maxiters = 100,
        kwarg...
)
    @unpack grid, imhomo, fun_para, ic, bc_para, Ny = prob
    _, nt, nx, ny = size(grid)
    Type = eltype(ic)
    nv = size(ic)[1] # nv is the number of solution vectors

    atol = abstol !== nothing ? abstol :
           real(oneunit(eltype(Type))) * (eps(real(one(eltype(Type)))))^(4 // 5)
    rtol = reltol !== nothing ? reltol : eps(real(one(eltype(Type))))^(4 // 5)

    flowS = complex(zeros(Type, nv, nt, nx, Ny))
    flowS0 = complex(zeros(Type, nv, nx, Ny))

    flow = zeros(Type, nv, nt, nx, ny)

    nonlinear_term = complex(zeros(Type, nv, nx, Ny))

    bcparaS = complex(similar(bc_para))
    icS = complex(similar(ic))

    ##* FFT for the forcing, bc and ic ------------------------------
    for i in 1:nt, n in 1:size(bc_para)[1]
        bcparaS[n, i, 1:Ny] = fft_expand(complex(bc_para[n, i, :]), Ny; atol = 1e-9)
    end
    for j in 1:nx, n in 1:nv
        icS[n, j, 1:Ny] = fft_expand(complex(ic[n, j, :]), Ny; atol = 1e-9)
    end

    ##* NSPPDE iteration ---------------------------------------------
    iterate && (println("NSPPDE iteration:"); showprogress = false)
    showprogress ? p = Progress(nt, "Solve the nonlinear PDE using spectral method") :
    nothing
    for i in 1:nt
        if i == 1
            flowS[:, i, :, :] = icS[:, :, 1:Ny]
            for n in 1:nv, j in 1:nx
                flow[n, i, j, :] = ifft_expand(flowS[n, i, j, :], ny)
            end
            showprogress && next!(p)
            continue
        else
            flowS0 = flowS[:, i - 1, :, :]
            flow[:, i, :, :] = flow[:, i - 1, :, :]
            for n in 1:maxiters
                imhomo(nonlinear_term, [flow, fun_para, Ny, i])
                for k in 1:Ny
                    M, B = jac(prob, flowS, nonlinear_term, bcparaS, i, k)
                    B = M \ B
                    for j in 1:nx
                        flowS[:, i, j, k] = B[((j - 1) * nv + 1):(j * nv)]
                    end
                end
                # IFFT for the solution
                for n in 1:nv, j in 1:nx
                    flow[n, i, j, :] = ifft_expand(flowS[n, i, j, :], ny)
                end
                if iterate
                    @printf "time = %e, iteration = %i, BW = %.3e\n" grid[1, i, 1, 1] n maximum(
                        abs.(flowS[:, i, :, :] - flowS0),
                    )
                else
                    break
                end
                isapprox(flowS0, flowS[:, i, :, :], atol = atol, rtol = rtol) && break
                flowS0 = flowS[:, i, :, :]
            end
            showprogress && next!(p)
        end
    end
    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function jac(prob, alg::PDEAlgorithm, flow, tn)
    @unpack fun, fun_para, imhomo, im_para, bc, bc_para, grid, ic = prob
    nx = size(grid)[3]
    T = eltype(ic)
    nv = size(ic)[1]

    M = zeros(T, nx * nv, nx * nv)
    B = zeros(T, nx * nv)

    Vxx, A, Γ = [zeros(T, nv, nv) for _ in 1:3]

    ΔT = grid[1, tn, 1] - grid[1, tn - 1, 1]

    if alg isa PSFDM
        Diff, x = chebyshevshift(nx, (grid[2, tn, 1], grid[2, tn, end]))
        Diff2 = Diff .^ 2
    elseif alg isa PFDM
        x = grid[2, tn, :]
        Diff, Diff2 = FDM_D(x)
    else
        error("Choose alg as 'PSFDM' or 'PFDM'. ")
    end

    for i in eachindex(x)
        ri = (1:nv) .+ ((i - 1) * nv)
        D = @view M[ri, ri]
        if fun_para isa NullParameter || ndims(fun_para) == 1
            fun(Vxx, A, Γ, D, fun_para)
        else
            fun(Vxx, A, Γ, D, fun_para[:, tn, i])
        end

        F = @view B[ri]
        if im_para isa NullParameter || isnothing(im_para) || ndims(im_para) == 1
            imhomo(F, im_para)
        else
            imhomo(F, im_para[:, tn, i])
        end

        if tn == 2
            C0 = 1
            Γtn0 = Γ / ΔT * flow[:, tn - 1, i]
        else
            C0 = 3 / 2
            Γtn0 = Γ / 2 / ΔT * (4 * flow[:, tn - 1, i] - flow[:, tn - 2, i])
        end
        M[ri, ri] += C0 .* Γ ./ ΔT
        B[ri] += Γtn0

        for j in eachindex(x)
            rj = (1:nv) .+ ((j - 1) * nv)
            M[ri, rj] += Diff2[i, j] * Vxx + Diff[i, j] * A
        end
    end

    bc(M, B, bc_para[:, tn])

    return M, B
end

function jac(prob, flow, imhomo, bc_para, tN, yN)
    @unpack fun, fun_para, bc, grid, ic = prob
    nx = size(grid)[3]
    T = eltype(ic)
    nv = size(ic)[1]

    M = zeros(T, nx * nv, nx * nv)
    B = zeros(T, nx * nv)

    Vxx, A, Γ, D = [zeros(T, nv, nv) for _ in 1:4]

    ΔT = grid[1, tN, 1, yN] - grid[1, tN - 1, 1, yN]

    Diff, x = chebyshevshift(nx, (grid[2, tN, 1, yN], grid[2, tN, end, yN]))
    grid[2, tN, :, yN] .= x

    for i in eachindex(x)
        ri = (1:nv) .+ ((i - 1) * nv)
        D = @view M[ri, ri]
        fun(Vxx, A, Γ, D, fun_para[:, tN, i, yN])

        B[ri] = imhomo[:, i, yN]

        if tN == 2
            C0 = 1
            Γtn0 = Γ / ΔT * flow[:, tN - 1, i, yN]
        else
            C0 = 3 / 2
            Γtn0 = Γ / 2 / ΔT * (4 * flow[:, tN - 1, i, yN] - flow[:, tN - 2, i, yN])
        end
        M[ri, ri] += C0 .* Γ ./ ΔT
        B[ri] += Γtn0
        for j in eachindex(x)
            rj = (1:nv) .+ ((j - 1) * nv)
            M[ri, rj] += Diff[i, j]^2 * Vxx + Diff[i, j] * A
        end
    end

    bc(M, B, bc_para[:, tN, yN])

    return M, B
end
