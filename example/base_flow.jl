using HydrodynamicStability, Plots
plotlyjs()

"""
    solve the base flow of the compressible blasius equations
"""

include("type.jl")
include("equations.jl");
# parameter setting
Ma = 3
Te = 250
Re = Inf
fs = FreeStream(Ma, Te, Re)

u₀ = [0, 0, 0.3779, 0, 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2, 0]
tspan = (0.0, 10.0)

# solve the compressible blasius equations in the coodinate η = 1 / sqrt(x * Re) ∫_0^y ρ(η) dη
bvp = BVProblem(similarity!, similarity_bc!, u₀, tspan, fs)
sol_B = solve(bvp, Shooting(), dy = 0.05);

## solve downstream flow
nx = 201
ny = 201
nz = 64
lx = 5
ly = 10

mesh = zeros(2, nx, ny)
flow = zeros(4, nx, ny)
ic = zeros(ny)

for i = 1:nx, j in eachindex(sol_B.y)
    mesh[1, i, j] = (i - 1) * lx / (nx - 1)
    mesh[2, i, j] = sol_B.y[j]
end

for j in eachindex(sol_B.y)
    ic[j] = sol_B.u[j][5]
end

M = ones(nx, ny)
A = ones(nx, ny)
D = zeros(nx, ny)
γ = ones(nx, ny)
bc_wall, bc_farfiled = [zeros(nx), zeros(nx)]

for i = 1:nx, j = 1:ny
    M[i, j] = 1 / Pr
    A[i, j] = sol_B.u[j][1]
    γ[i, j] = -2 * mesh[1, i, j] * sol_B.u[j][2]
end

h = 0.3
d = 0.3
x_0 = 0.955 # The maximum is located at x = (x_0 + √(x_0^2 + 2d^2)) / 2, when d = 0.3 x_0 = 0.955, x_max = 1 
F(x) = x * exp(-(x - x_0)^2 / d^2)
T_B0 = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2

for i = 1:nx
    bc_wall[i] = T_B0
    bc_farfiled[i] = 1
end

f_para = [M, A, D, γ]
bc_para = [bc_wall, bc_farfiled]

prob = HeatProblem(mesh, ic, f_para, bc_para)
sol = solve(prob, FDM())
