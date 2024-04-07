using HydrodynamicStability, Test

"""
    solve the Heat equation for the test of PDEProblem
"""

nt = 201
nx = 201
lt = 1
lx = 1

mesh = zeros(2, nt, nx)

for i = 1:nt, j = 1:nx
    mesh[1, i, j] = (i - 1) * lt / (nt - 1)
    mesh[2, i, j] = (j - 1) * lx / (nx - 1)
end

function heat_fun!(Vxx, A, Γ, D, F, p)
    λ, = p
    A[1, 1] = 1
    A[2, 2] = 1
    Γ[2, 1] = -λ
    D[1, 2] = -1
end

function heat_bc!(M0, Mend, u, p)
    x₀, x∞ = p
    M0[2, :] .= 0
    M0[2, 1] = 1
    Mend[2, :] .= 0
    Mend[2, end-1] = 1

    u[2] = x₀
    u[end] = x∞
end

ic = Array{Array{Float64}}(undef, nx)
for i = 1:nx
    ic[i] =
        [sin(2π * mesh[2, 1, i])^2, 4π * cos(2π * mesh[2, 1, i]) * sin(2π * mesh[2, 1, i])]
end

f_para = Array{Array{Float64}}(undef, nt, nx)
bc_para = Array{Array{Float64}}(undef, nt)

for i = 1:nt
    for j = 1:nx
        f_para[i, j] = [1]
    end
    bc_para[i] = [0, 0]
end

prob = HeatProblem(mesh, heat_fun!, heat_bc!, ic, f_para, bc_para)
sol = solve(prob, FDM())

@test sol.flow[1, 1, :] == sin.(2π * mesh[2, 1, :]) .^ 2
@test isapprox(maximum(sol.flow[1, end, :]), 0.0, atol = 1e-4)

# TODO find the reason why the result is oscillating, Runga's phenomenon ?
