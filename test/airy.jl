using HydrodynamicStability, Test, SpecialFunctions

"""
    Airy_fun!(D, F, p, t)
    
solve the Airy equation
```math
u'' - t u = 0,
```
with the boundary condition
```math
u(0) = Ai(0)
u(10) = Ai(10),
```
where ``Ai(x)`` is the Airy function of the first kind.

"""
function AiryFun!(A, D, F, t, p)
    A[1, 1] = 1
    A[2, 2] = 1

    D[1, 2] = -1
    D[2, 1] = -t
end

function AiryBC!(M, u)
    M[2, :] .= 0
    M[2, 1] = 1
    M[end, :] .= 0
    M[end, end - 1] = 1

    u[2] = airyai(-10)
    u[end] = airyai(10)
end

yspan = range(-10, 10, 2001)
u₀ = [0.0, 0.0]
prob = BVProblem(AiryFun!, AiryBC!, u₀, yspan)
ret = solve(prob, FDM(); ny = 2001)
sol = flatten_vector(ret.u)

yrange = (-10, 10)
probc = BVProblem(AiryFun!, AiryBC!, u₀, yrange)
retc = solve(probc, SFDM(); ny = 64)
solc = flatten_vector(retc.u)

@test isapprox(sqrt(sum(abs2, sol[2, :] - airyaiprime.(ret.y))), 0.0, atol = 1e-3)
@test isapprox(solc[2, :], airyaiprime.(retc.y))

##* Nonlinear Equation ----------------------------------------------

function FunProblrm1!(A, D, F, t, p)
    A[1, 1] = 1
    A[2, 2] = 1

    D[1, 2] = -1
    D[2, 1] = -1
end

function NProblem1!(N, u, t, p)
    N[2] = -u[1]^2
end

function DNProblem1!(DN, u, t, p)
    DN[2, 1] = -2 * u[1]
end

funs = (FunProblrm1!, NProblem1!, DNProblem1!)

function BCProblem1!(M, u)
    M[2, :] .= 0
    M[2, 1] = 1
    M[end, :] .= 0
    M[end, end - 1] = 1

    u[2] = 1
    u[end] = 4
end

function NBCProblem1!(N)
    N[2] = 0
    N[end] = 0
end

function DNBCProblem1!(DN)
    DN[2, :] .= 0
    DN[end, :] .= 0
end

bcs = (BCProblem1!, NBCProblem1!, DNBCProblem1!)

yspan = (0, 1)
ny = 64
u₀ = ones(2 * ny)

prob = BVProblem(funs, bcs, u₀, yspan)
ret = solve(prob, NSFDM(); ny = ny)
sol = flatten_vector(ret.u)

D, tspan = chebyshevshift(ny, yspan)
fun1 = D * sol[1, :] - sol[2, :]
fun2 = D * sol[2, :] - (sol[1, :] - sol[1, :] .^ 2)
@test isapprox(sqrt(sum(abs2, fun1)), 0, atol = 1e-9)
@test isapprox(sqrt(sum(abs2, fun2)), 0, atol = 1e-9)

##* Blasius Equation ------------------------------------------------

# function Blasius!(A, D, F, t, p)
#     A[1, 1] = 1
#     A[2, 2] = 1
#     A[3, 3] = 1

#     D[1, 2] = -1
#     D[2, 3] = -1
# end

# function NBlasius!(N, u, t, p)
#     N[3] = -1 / 2 * u[1] * u[3]
# end

# function DNBlasius!(DN, u, t, p)
#     DN[3, 1] = -1 / 2 * u[3]
#     DN[3, 3] = -1 / 2 * u[1]
# end

# funs = (Blasius!, NBlasius!, DNBlasius!)

# function BlasiusBC!(M, u)
#     M[1, :] .= 0
#     M[1, 1] = 1
#     M[3, :] .= 0
#     M[3, 2] = 1
#     M[end, :] .= 0
#     M[end, end - 1] = 1

#     u[1] = 0
#     u[3] = 0
#     u[end] = 1
# end

# function NBlasiusBC!(N)
#     N[1] = 0
#     N[3] = 0
#     N[end] = 0
# end

# function DNBlasiusBC!(DN)
#     DN[1, :] .= 0
#     DN[3, :] .= 0
#     DN[end, :] .= 0
# end

# bcs = (BlasiusBC!, NBlasiusBC!, DNBlasiusBC!)

# yspan = (0, 15)
# ny = 64
# # u₀ = 1 / 3 * ones(3 * ny)
# u₀ = u0

# prob = BVProblem(funs, bcs, u₀, yspan)
# ret = solve(prob, NSFDM(); ny = ny, showiters = true)
# sol = flatten_vector(ret.u)

# fun1 = D * sol[1, :] - sol[2, :]
# fun2 = D * sol[2, :] - sol[3, :]
# fun3 = D * sol[3, :] + 1 / 2 * sol[1, :] .* sol[3, :]
# @test isapprox(sqrt(sum(abs2, fun1)), 0, atol = 1e-9)
# @test isapprox(sqrt(sum(abs2, fun2)), 0, atol = 1e-9)
# @test isapprox(sqrt(sum(abs2, fun3)), 0, atol = 1e-9)
