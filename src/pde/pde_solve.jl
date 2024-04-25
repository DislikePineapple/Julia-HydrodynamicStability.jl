"""
    solve(prob::PDEProblem, alg::FDM, args...; kwarg...)

Solve the PDE problem using the backward finite difference method with fixed step-length grid.
```math
\frac{∂^2u}{∂x^2} - λ\frac{∂u}{∂t} = 0 
```
"""

# using FiniteDiff to solve the parabolic PDE problem
abstract type PPDEAlgorithm <: AbstractAlgorithm end
abstract type NPPDEAlgorithm <: PPDEAlgorithm end
struct NSFDM <: NPPDEAlgorithm end # nonlinear spectral finite difference method

function solve(prob::HeatProblem, alg::FDM, args...; kwarg...)
    @unpack grid, ic = prob
    _, nt, nx = size(grid)
    Type = eltype(ic)
    nv = size(ic)[1]

    flow = zeros(Type, nv, nt, nx)

    for j = 1:nt
        if j == 1
            flow[:, j, :] = ic
            continue
        else
            M, B = jac(prob, flow, j)
            B = M \ B
            for k = 1:nx
                flow[:, j, k] = B[(k-1)*nv+1:k*nv]
            end
        end
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function solve(
    prob::NSHeatProblem,
    alg::NSFDM,
    args...;
    abstol = nothing,
    reltol = nothing,
    showiters = false,
    maxiters = 100,
    kwarg...,
)
    @unpack grid, imhomo, imhomo_para, ic, bc_para, Ny = prob
    _, nt, nx, ny = size(grid)
    Type = eltype(ic)
    nv = size(ic)[1] # nv is the number of solution vectors

    atol =
        abstol !== nothing ? abstol :
        real(oneunit(eltype(Type))) * (eps(real(one(eltype(Type)))))^(4 // 5)
    rtol = reltol !== nothing ? reltol : eps(real(one(eltype(Type))))^(4 // 5)

    flowS = complex(zeros(Type, nv, nt, nx, Ny))
    flowS0 = complex(zeros(Type, nv, nx, Ny))

    bcparaS = complex(similar(bc_para))
    icS = complex(similar(ic))

    ##* FFT for the forcing, bc and ic ------------------------------
    for i = 1:nt, n = 1:size(bc_para)[1]
        bcparaS[n, i, 1:Ny] = fft_expand(complex(bc_para[n, i, :]), Ny; atol = 1e-9)
    end
    for j = 1:nx, n = 1:nv
        icS[n, j, 1:Ny] = fft_expand(complex(ic[n, j, :]), Ny; atol = 1e-9)
    end

    ##* NPPDE iteration ---------------------------------------------
    showiters && println("NPPDE iteration:")
    tspan = grid[1, :, 1, 1]
    for i = 1:nt
        if i == 1
            flowS[:, i, :, :] = icS[:, :, 1:Ny]
            continue
        else
            flowS0[:, :, :] = flowS[:, i-1, :, :]
            for n = 1:maxiters
                imhomoS = imparaS(imhomo, flowS0, imhomo_para, ny)
                for k = 1:Ny
                    M, B = jac(prob, flowS, imhomoS, bcparaS, i, k)
                    B = M \ B
                    for j = 1:nx
                        flowS[:, i, j, k] = B[(j-1)*nv+1:j*nv]
                    end
                end
                showiters &&
                    @printf "time = %e, iteration = %i, BW = %.3e\n" tspan[i] n maximum(
                        abs.(flowS[:, i, :, :] - flowS0),
                    )
                isapprox(flowS0, flowS[:, i, :, :], atol = atol, rtol = rtol) && break
                flowS0 = flowS[:, i, :, :]
            end
        end
    end

    ##* IFFT for the solution ---------------------------------------
    flow = zeros(Type, nv, nt, nx, ny)
    for n = 1:nv, i = 1:nx, j = 1:nx
        flow[n, i, j, :] = ifft_expand(flowS[n, i, j, :], ny)
    end

    HeatSolution(grid, flow, prob, alg, nothing, nothing)
end

function jac(prob, flow, tn)
    @unpack fun, fun_para, imhomo, im_para, bc, bc_para, grid, ic = prob
    nx = size(grid)[3]
    T = eltype(ic)
    nv = size(ic)[1]

    M = zeros(T, nx * nv, nx * nv)
    B = zeros(T, nx * nv)

    Vxx, A, Γ, D = [zeros(T, nv, nv) for _ = 1:4]
    F = zeros(T, nv)

    x = grid[2, tn, :]
    ΔT = grid[1, tn, 1] - grid[1, tn-1, 1]

    tn == 2 ? C0 = 1 : C0 = 3 / 2

    for i in eachindex(x)
        position, dx = fun_position(x, i)
        for j in eachindex(x)
            if fun_para isa NullParameter || ndims(fun_para) == 1
                fun(Vxx, A, Γ, D, fun_para)
            else
                fun(Vxx, A, Γ, D, fun_para[:, tn, j])
            end
            if im_para isa NullParameter || isnothing(im_para) || ndims(im_para) == 1
                imhomo(F, im_para)
            else
                imhomo(F, im_para[:, tn, j])
            end

            if i == j
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 3] .* Vxx ./ dx^2 +
                    C1[position, 3] .* A ./ dx +
                    D +
                    C0 .* Γ ./ ΔT
                if tn == 2
                    B[(i-1)*nv+1:i*nv] = Γ * flow[:, tn-1, i] ./ ΔT + F
                else
                    B[(i-1)*nv+1:i*nv] =
                        Γ * (4 .* flow[:, tn-1, i] - flow[:, tn-2, i]) ./ 2 ./ ΔT + F
                end
            elseif j == i - 2
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 1] .* Vxx ./ dx^2 + C1[position, 1] .* A ./ dx
            elseif j == i - 1
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 2] .* Vxx ./ dx^2 + C1[position, 2] .* A ./ dx
            elseif j == i + 1
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 4] .* Vxx ./ dx^2 + C1[position, 4] .* A ./ dx
            elseif j == i + 2
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 5] .* Vxx ./ dx^2 + C1[position, 5] .* A ./ dx
            end
        end
    end

    M0 = M[1:nv, :]
    Mend = M[end-nv+1:end, :]

    bc(M0, Mend, B, bc_para[:, tn])

    M[1:nv, :] = M0
    M[end-nv+1:end, :] = Mend

    return M, B
end

function jac(prob, flow, imhomo, bc_para, tN, yN)
    @unpack fun, fun_para, bc, grid, ic = prob
    nx = size(grid)[3]
    T = eltype(ic)
    nv = size(ic)[1]

    M = zeros(T, nx * nv, nx * nv)
    B = zeros(T, nx * nv)

    Vxx, A, Γ, D = [zeros(T, nv, nv) for _ = 1:4]
    F = zeros(T, nv)

    x = grid[2, tN, :, yN]
    ΔT = grid[1, tN, 1, yN] - grid[1, tN-1, 1, yN]

    tN == 2 ? C0 = 1 : C0 = 3 / 2

    for i in eachindex(x)
        position, dx = fun_position(x, i)
        for j in eachindex(x)
            fun(Vxx, A, Γ, D, fun_para[:, tN, j, yN])
            F = imhomo[:, j, yN]
            if i == j
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 3] .* Vxx ./ dx^2 +
                    C1[position, 3] .* A ./ dx +
                    D +
                    C0 .* Γ ./ ΔT
                if tN == 2
                    B[(i-1)*nv+1:i*nv] = Γ * flow[:, tN-1, i, yN] ./ ΔT + F
                else
                    B[(i-1)*nv+1:i*nv] =
                        Γ * (4 .* flow[:, tN-1, i, yN] - flow[:, tN-2, i, yN]) ./ 2 ./ ΔT +
                        F
                end
            elseif j == i - 2
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 1] .* Vxx ./ dx^2 + C1[position, 1] .* A ./ dx
            elseif j == i - 1
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 2] .* Vxx ./ dx^2 + C1[position, 2] .* A ./ dx
            elseif j == i + 1
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 4] .* Vxx ./ dx^2 + C1[position, 4] .* A ./ dx
            elseif j == i + 2
                M[(i-1)*nv+1:i*nv, (j-1)*nv+1:j*nv] +=
                    C2[position, 5] .* Vxx ./ dx^2 + C1[position, 5] .* A ./ dx
            end
        end
    end

    M0 = M[1:nv, :]
    Mend = M[end-nv+1:end, :]

    bc(M0, Mend, B, bc_para[:, tN, yN])

    M[1:nv, :] = M0
    M[end-nv+1:end, :] = Mend

    return M, B
end

function imparaS(imhomo!, flowS, coe, ny)
    nv, nx, Ny = size(flowS)
    flow = zeros(nv, nx, ny)
    #Fourier inverse transform
    for n = 1:nv, j = 1:nx
        flow[n, j, :] = ifft_expand(flowS[n, j, :], ny)
        # selfdefined function `ifft_expand` in "utils.jl"
    end

    imhomo = zeros(nv)
    imhomoS = zeros(ComplexF64, nv, nx, ny)
    for j = 1:nx, k = 1:ny
        imhomo!(imhomo, (coe, flow[:, j, k]))
        imhomoS[:, j, k] = imhomo
    end

    for n = 1:nv, j = 1:nx
        imhomoS[n, j, 1:Ny] = fft_expand(imhomoS[n, j, :], Ny; atol = 1e-3)
    end
    return imhomoS[:, :, 1:Ny]
end
