using HydrodynamicStability, Test, SpecialFunctions

"""
    Airy_fun!(M, B, p, t)
    solve the Airy equation
    u'' - t * u = 0,
    which can be write as the first order ODE system
    u' = v
    v' = t * u,
    with the boundary condition
    u(0) = Ai(0)
    u(10) = Ai(10),
    where Ai(x) is the Airy function of the first kind.
"""

function Airy_fun!(A, D, F, p, t)
    A[1, 1] = 1
    A[2, 2] = 1

    D[1, 2] = -1
    D[2, 1] = -t

    F[1] = 0
    F[2] = 0
end

function Airy_bc!(M0, Mend, u)
    M0[2, :] .= 0
    M0[2, 1] = 1
    Mend[2, :] .= 0
    Mend[2, end-1] = 1

    u[2] = airyai(0)
    u[end] = airyai(10)
end

y = 0:0.01:10
u₀ = [0.0, 0.0]
prob = BVProblem(Airy_fun!, Airy_bc!, u₀, y)
sol = solve(prob, FDM())

@test isapprox(sol.u[1][2], airyaiprime(0), atol = 1e-5)
