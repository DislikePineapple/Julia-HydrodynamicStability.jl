function solve(prob::ODEProblem, alg::RK4, args...; dy, kwargs...)
    @unpack f, u0, yspan, p = prob
    y = collect(yspan[1]:dy:yspan[2])

    u = Array{typeof(u0)}(undef, length(y))

    u[1] = u0

    s1 = similar(u0)
    s2 = similar(u0)
    s3 = similar(u0)
    s4 = similar(u0)

    for i = 1:length(y)-1
        h = y[i+1] - y[i]
        f(s1, u[i], p, y[i])
        f(s2, u[i] .+ h / 2 * s1, p, y[i] + h / 2)
        f(s3, u[i] .+ h / 2 * s2, p, y[i] + h / 2)
        f(s4, u[i] .+ h * s3, p, y[i] + h)
        u[i+1] = u[i] .+ h / 6 * (s1 .+ 2 * s2 .+ 2 * s3 .+ s4)
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

function solve(prob::BVProblem, alg::FDM, args...; kwarg...)
    @unpack yspan, u0 = prob
    T = eltype(u0)
    O = length(u0)

    y = collect(yspan)
    u = Array{Array{T}}(undef, length(y))

    M, B = jac(prob)
    B = M \ B

    for i = 1:length(y)
        u[i] = B[(i-1)*O+1:i*O]
    end
    ODESolution(u, y)
end

function jac(prob::BVProblem)
    @unpack f, bc, p, yspan, u0 = prob
    T = eltype(u0)
    O = length(u0)

    y = collect(yspan)

    ny = length(y)

    M = zeros(T, ny * O, ny * O)
    B = zeros(T, ny * O)

    A, D = [zeros(T, O, O) for _ = 1:2]
    F = zeros(T, O)

    for i in eachindex(y)

        position, dy = fun_position(y, i)

        for j in eachindex(y)
            if p isa NullParameter || p[1] isa Number
                f(A, D, F, p, y[i])
            else
                f(A, D, F, p[i], y[i])
            end
            if i == j
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 3] .* A ./ dy + D
                B[(i-1)*O+1:i*O] = F
            elseif j == i - 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 1] .* A ./ dy
            elseif j == i - 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 2] .* A ./ dy
            elseif j == i + 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 4] .* A ./ dy
            elseif j == i + 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 5] .* A ./ dy
            end
        end
    end

    M0 = M[1:O, :]
    Mend = M[end-O+1:end, :]

    bc(M0, Mend, B)

    M[1:O, :] = M0
    M[end-O+1:end, :] = Mend

    return M, B
end

function fun_position(y, i)
    if i == 1
        position = 1
        dy = y[i+1] - y[i]
    elseif i == 2
        position = 2
        dy = (y[i+1] - y[i-1]) / 2
    elseif i == length(y) - 1
        position = 4
        dy = (y[i+1] - y[i-1]) / 2
    elseif i == length(y)
        position = 5
        dy = y[i] - y[i-1]
    else
        position = 3
        dy = (y[i+1] - y[i-1]) / 2
    end
    return position, dy
end
