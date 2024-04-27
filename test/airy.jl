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

function Airy_fun!(D, F, p, t)
    D[1, 2] = -1
    D[2, 1] = -t
end

function Airy_bc!(M, u)
    M[2, :] .= 0
    M[2, 1] = 1
    M[end, :] .= 0
    M[end, end - 1] = 1

    u[2] = airyai(0)
    u[end] = airyai(10)
end

yspan = (0, 10)
u₀ = [0.0, 0.0]
prob = BVProblem(Airy_fun!, Airy_bc!, u₀, yspan)

sol = solve(prob, FDM(); ny = 1001)
csol = solve(prob, CFDM(); ny = 32)

@test isapprox(sol.u[2, 1], airyaiprime(0), atol = 1e-6)
@test isapprox(csol.u[2, 1], airyaiprime(0))
