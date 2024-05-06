function solve(prob::ODEProblem, alg::RK4, args...; kwargs...)
    @unpack f, u0, yspan, p = prob

    y = collect(yspan)
    u = Array{typeof(u0)}(undef, length(y))
    s1, s2, s3, s4 = [similar(u0) for _ in 1:4]

    u[1] = u0
    for i in 1:(length(y) - 1)
        h = y[i + 1] - y[i]
        f(s1, u[i], y[i], p)
        f(s2, u[i] .+ h / 2 * s1, y[i] + h / 2, p)
        f(s3, u[i] .+ h / 2 * s2, y[i] + h / 2, p)
        f(s4, u[i] .+ h * s3, y[i] + h, p)
        u[i + 1] = u[i] .+ h / 6 * (s1 .+ 2 * s2 .+ 2 * s3 .+ s4)
    end

    return ODESolution(u, y)
end

function solve(prob::BVProblem, alg::Shooting, args...; kwargs...)
    @unpack f, bc, yspan, u0, p = prob
    @unpack ivp, iter = alg

    y = collect(yspan)
    bc!(residual, u) = bc(residual, u, y, p)

    function f_non(u0, p)
        residual = similar(u0)
        sol = solve(ODEProblem(f, u0, y, p), ivp)
        bc!(residual, sol.u)
        return residual
    end
    ic = solve(NonlinearProblem(f_non, u0; p = p); kwargs...)

    return solve(ODEProblem(f, ic.t, y, p), ivp)
end

function solve(prob::BVProblem, alg::FDM, args...; kwarg...)
    @unpack u0, yspan = prob
    nv = length(u0)

    M, B = jac(prob)
    B = M \ B

    u = nest_vector(B, nv)
    return ODESolution(u, yspan)
end

function jac(prob::BVProblem)
    @unpack f, bc, p, u0, yspan = prob
    T = eltype(u0) # type of solution vector
    nv = length(u0) # dims of solution vector
    ny = length(yspan) # number of grid points

    M = zeros(T, ny * nv, ny * nv) # Jacobian matrix
    B = zeros(T, ny * nv) # inhomogeneous vector

    Diff = FDM_D(yspan)[1] # finite difference matrix

    # Assembly matrix and vector
    A = zeros(T, nv, nv)
    for i in eachindex(yspan)
        ri = (1:nv) .+ ((i - 1) * nv)
        D = @view M[ri, ri]
        F = @view B[ri]
        if p isa NullParameter || p[1] isa Number
            f(A, D, F, yspan[i], p)
        else
            f(A, D, F, yspan[i], p[i])
        end
        for j in eachindex(yspan)
            rj = (1:nv) .+ ((j - 1) * nv)
            M[ri, rj] += Diff[i, j] * A
        end
    end

    bc(M, B) # apply boundary conditions
    return M, B
end

function solve(prob::BVProblem, alg::SFDM, args...; ny, kwarg...)
    @unpack yspan, u0 = prob
    nv = length(u0) # dims of solution vector

    D, y = chebyshevshift(ny, yspan) # Chebyshev grid, and rescale to yspan

    M, B = jac(prob, alg, D, y) # Jacobian matrix and inhomogeneous vector
    B = M \ B # solve the linear system

    u = nest_vector(B, nv) # reshape the solution vector
    return ODESolution(u, y)
end

function solve(prob::BVProblem, alg::NSFDM, args...; ny, kwargs...)
    @unpack yspan, u0 = prob
    nv = Int(length(u0) / ny) # dims of solution vector

    # Chebtshev grid, and rescale to the interval tspan
    D, y = chebyshevshift(ny, yspan)
    F, DF = jac(prob, alg, D, y)

    # return F, DF

    prob = NonlinearProblem((F, DF), u0)
    sol = solve(prob, Newton(); kwargs...)

    u = nest_vector(sol.t, nv)
    return ODESolution(u, y)
end

function jac(prob::BVProblem, alg::SFDM, Diff::AbstractArray, y::AbstractArray)
    @unpack f, bc, p, u0 = prob
    T = eltype(u0) # type of solution vector
    ny = length(y) # number of grid points
    nv = length(u0) # dims of solution vector

    M = zeros(T, ny * nv, ny * nv)
    B = zeros(T, ny * nv)

    A = zeros(T, nv, nv)
    for i in eachindex(y)
        ri = (1:nv) .+ ((i - 1) * nv)
        D = @view M[ri, ri]
        F = @view B[ri]
        if p isa NullParameter || p[1] isa Number
            f(A, D, F, y[i], p)
        else
            f(A, D, F, y[i], p[i])
        end
        for j in eachindex(y)
            rj = (1:nv) .+ ((j - 1) * nv)
            M[ri, rj] += Diff[i, j] * A
        end
    end

    bc(M, B)
    return M, B
end

function jac(prob::BVProblem, alg::NSFDM, Diff::AbstractArray, y::AbstractArray)
    @unpack f, bc, p, u0 = prob
    T = eltype(u0) # type of solution vector
    ny = length(y)
    nv = Int(length(u0) / ny) # dims of solution vector

    M = zeros(T, ny * nv, ny * nv)
    B = zeros(T, ny * nv)

    A = zeros(T, nv, nv)

    fun = prob.f[1]
    Nfun = prob.f[2]
    DNfun = prob.f[3]

    bc = prob.bc[1]
    Nbc = prob.bc[2]
    DNbc = prob.bc[3]

    for i in eachindex(y)
        ri = (1:nv) .+ ((i - 1) * nv)
        D = @view M[ri, ri]
        F = @view B[ri]
        if p isa NullParameter || p[1] isa Number
            fun(A, D, F, y[i], p)
        else
            fun(A, D, F, y[i], p[i])
        end
        for j in eachindex(y)
            rj = (1:nv) .+ ((j - 1) * nv)
            M[ri, rj] += Diff[i, j] * A
        end
    end

    bc(M, B)

    n = zeros(T, ny * nv)
    dn = zeros(T, ny * nv, ny * nv)

    function F(u, p)
        for i in eachindex(y)
            ri = (1:nv) .+ ((i - 1) * nv)
            N = @view n[ri]
            if p isa NullParameter || p[1] isa Number
                Nfun(N, u[ri], y[i], p)
            else
                Nfun(N, u[ri], y[i], p[i])
            end
        end
        Nbc(n)
        return M * u - n - B
    end

    function DF(du, u, p)
        for i in eachindex(y)
            ri = (1:nv) .+ ((i - 1) * nv)
            DN = @view dn[ri, ri]
            DNfun(DN, u[ri], y[i], p)
        end
        DNbc(dn)
        return M - dn
    end

    return F, DF
end
