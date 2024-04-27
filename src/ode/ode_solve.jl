function solve(prob::ODEProblem, alg::RK4, args...; dy, kwargs...)
    @unpack f, u0, yspan, p = prob
    y = collect(yspan[1]:dy:yspan[2])

    u = Array{typeof(u0)}(undef, length(y))

    u[1] = u0

    s1, s2, s3, s4 = [similar(u0) for _ in 1:4]

    for i in 1:(length(y) - 1)
        h = y[i + 1] - y[i]
        f(s1, u[i], p, y[i])
        f(s2, u[i] .+ h / 2 * s1, p, y[i] + h / 2)
        f(s3, u[i] .+ h / 2 * s2, p, y[i] + h / 2)
        f(s4, u[i] .+ h * s3, p, y[i] + h)
        u[i + 1] = u[i] .+ h / 6 * (s1 .+ 2 * s2 .+ 2 * s3 .+ s4)
    end

    ODESolution(u, y)
end

function solve(prob::BVProblem, alg::Shooting, args...; dy, kwargs...)
    @unpack f, bc, yspan, u0, p = prob
    @unpack ivp, iter = alg

    y = collect(yspan[1]:dy:yspan[2])
    bc!(residual, u) = bc(residual, u, p, y)

    function f_non(u0, p)
        residual = similar(u0)
        sol = solve(ODEProblem(f, u0, yspan, p), ivp, dy = dy)
        bc!(residual, sol.u)
        return residual
    end

    ic = solve(NonlinearProblem(f_non, u0, p); kwargs...)
    solve(ODEProblem(f, ic.t, yspan, p), ivp, dy = dy)
end

function solve(prob::BVProblem, alg::CFDM, args...; ny, kwarg...)
    @unpack yspan, u0 = prob
    T = eltype(u0) # type of solution vector
    nv = length(u0) # dims of solution vector

    # 1D Chebyshev grid, and rescale to the interval yspan
    D, y = chebyshev(ny)
    y = yspan[1] .+ (yspan[2] - yspan[1]) .* (y .+ 1) ./ 2
    D = D ./ ((yspan[2] - yspan[1]) / 2)

    M, B = jac(prob, D, y)
    B = M \ B

    u = zeros(T, nv, ny)
    for i in 1:ny
        u[:, i] = B[((i - 1) * nv + 1):(i * nv)]
    end
    ODESolution(u, y)
end

function solve(prob::BVProblem, alg::FDM, args...; ny, kwarg...)
    @unpack yspan, u0 = prob
    T = eltype(u0)
    nv = length(u0)

    y = range(yspan[1], yspan[2], length = ny)

    M, B = jac(prob, y)
    B = M \ B

    u = zeros(T, nv, ny)
    for i in 1:ny
        u[:, i] = B[((i - 1) * nv + 1):(i * nv)]
    end
    ODESolution(u, y)
end

function jac(prob::BVProblem, Diff::AbstractArray, y::AbstractArray)
    @unpack f, bc, p, u0 = prob
    T = eltype(u0) # type of solution vector
    nv = length(u0) # dims of solution vector
    ny = length(y)

    M = zeros(T, ny * nv, ny * nv)
    B = zeros(T, length(y) * nv)

    for i in eachindex(y)
        D = @view M[(1 + (i - 1) * nv):(i * nv), (1 + (i - 1) * nv):(i * nv)]
        F = @view B[(1 + (i - 1) * nv):(i * nv)]
        if p isa NullParameter || p[1] isa Number
            f(D, F, p, y[i])
        else
            f(D, F, p[i], y[i])
        end
    end
    M += kron(Diff, I(nv))
    bc(M, B)

    return M, B
end

function jac(prob::BVProblem, y::AbstractArray)
    @unpack f, bc, p, u0 = prob
    T = eltype(u0) # type of solution vector
    nv = length(u0) # dims of solution vector
    ny = length(y)

    M = zeros(T, ny * nv, ny * nv)
    B = zeros(T, length(y) * nv)

    for i in eachindex(y)
        D = @view M[(1 + (i - 1) * nv):(i * nv), (1 + (i - 1) * nv):(i * nv)]
        F = @view B[(1 + (i - 1) * nv):(i * nv)]
        if p isa NullParameter || p[1] isa Number
            f(D, F, p, y[i])
        else
            f(D, F, p[i], y[i])
        end
    end

    Diff = FDM_D(y)[1]
    M += kron(Diff, I(nv))
    bc(M, B)

    return M, B
end

function FDM_D(y)
    D, D2 = [zeros(length(y), length(y)) for _ in 1:2]
    for i in eachindex(y)
        if i == 1
            position = 1
            dy = y[i + 1] - y[i]
        elseif i == 2
            position = 2
            dy = (y[i + 1] - y[i - 1]) / 2
        elseif i == length(y) - 1
            position = 4
            dy = (y[i + 1] - y[i - 1]) / 2
        elseif i == length(y)
            position = 5
            dy = y[i] - y[i - 1]
        else
            position = 3
            dy = (y[i + 1] - y[i - 1]) / 2
        end
        for j in eachindex(y)
            if p isa NullParameter || p[1] isa Number
                f(A, D, F, p, y[i])
            else
                f(A, D, F, p[i], y[i])
            end
            if i == j
                D[i, j] = C1[position, 3] / dy
                D2[i, j] = C2[position, 3] / dy^2
            elseif j == i - 2
                D[i, j] = C1[position, 1] / dy
                D2[i, j] = C2[position, 1] / dy^2
            elseif j == i - 1
                D[i, j] = C1[position, 2] / dy
                D2[i, j] = C2[position, 2] / dy^2
            elseif j == i + 1
                D[i, j] = C1[position, 4] / dy
                D2[i, j] = C2[position, 4] / dy^2
            elseif j == i + 2
                D[i, j] = C1[position, 5] / dy
                D2[i, j] = C2[position, 5] / dy^2
            end
        end
    end
    return D, D2
end
