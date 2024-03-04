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
    @unpack f, bc, yspan, u0, p = prob
    T = eltype(u0)
    O = length(u0)

    y = collect(yspan)

    ny = length(y)
    u = Array{Array{T}}(undef, ny)

    M = zeros(T, ny * O, ny * O)
    B = zeros(T, ny * O)

    jac!(M, B, f, bc, p, y, O, T)
    B = M \ B
    for i = 1:ny
        u[i] = B[(i-1)*O+1:i*O]
    end
    ODESolution(u, y)
end

function jac!(M, B, f!, bc!, p, y, O, T)
    M1 = zeros(T, O, O)
    B1 = zeros(T, O)

    for i in eachindex(y)

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

        for j in eachindex(y)
            if i == j
                if p isa NullParameter || p[1] isa Number
                    f!(M1, B1, p, y[i])
                else
                    f!(M1, B1, p[i], y[i])
                end
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 3] * I(O) ./ dy - M1
                B[(i-1)*O+1:i*O] = B1
            elseif j == i - 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 1] * I(O) ./ dy
            elseif j == i - 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 2] * I(O) ./ dy
            elseif j == i + 1
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 4] * I(O) ./ dy
            elseif j == i + 2
                M[(i-1)*O+1:i*O, (j-1)*O+1:j*O] += C1[position, 5] * I(O) ./ dy
            end
        end
    end

    M0 = M[1:O, :]
    Mend = M[end-O+1:end, :]

    bc!(M0, Mend, B)

    M[1:O, :] = M0
    M[end-O+1:end, :] = Mend
end
