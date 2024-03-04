using HydrodynamicStability, Test

"""
    solve the Heat equation for the test of PDEProblem
"""

nt = 200
nx = 200
lt = 1
lx = 1

mesh = zeros(2, nt, nx)

for i = 1:nt, j = 1:nx
    mesh[1, i, j] = (i - 1) * lt / (nt - 1)
    mesh[2, i, j] = (j - 1) * lx / (nx - 1)
end

ic = sin.(2π * mesh[2, 1, :]) .^ 2
f_para = [ones(nt, nx), zeros(nt, nx), zeros(nt, nx), -ones(nt, nx)]
bc_para = [zeros(nt), zeros(nt)]

prob = HeatProblem(mesh, ic, f_para, bc_para)
sol = solve(prob, CrankNicolson())

@test sol.flow[1, 1, :] == ic
@test isapprox(maximum(sol.flow[1, end, :]), 0.0, atol = 1e-4)

# TODO find the reason why the result is oscillating, Runga's phenomenon ?
