using HydrodynamicStability, Test

"""
    solve the Heat equation for the test of PDEProblem
"""
nt = 101
nx = 101
ny = 101
lt = 1
lx = 1
lys = -π
lye = π

Ny = 4

xspan = range(0, lx, nx)
tspan = range(0, lt, nt)
mesh = zeros(2, nt, nx)

for i in 1:nt, j in 1:nx
    mesh[1, i, j] = (i - 1) * lt / (nt - 1)
    mesh[2, i, j] = (j - 1) * lx / (nx - 1)
end

function heat_fun!(Vxx, A, Γ, D, p)
    λ, = p
    A[1, 1] = 1
    A[2, 2] = 1
    Γ[2, 1] = -λ
    D[1, 2] = -1
end
function imhomo!(F, p) end
function heat_bc!(M, u, p)
    x₀, x∞ = p
    M[2, :] .= 0
    M[2, 1] = 1
    M[end, :] .= 0
    M[end, end - 1] = 1

    u[2] = x₀
    u[end] = x∞
end
funs = (heat_fun!, imhomo!, heat_bc!)

ic = zeros(2, nx)
for i in 1:nx
    ic[1, i] = sin(2π * xspan[i])^2
    ic[2, i] = 4π * cos(2π * xspan[i]) * sin(2π * xspan[i])
end
fp = ones(1, nt, nx)
bcp = zeros(2, nt)
paras = (fp, nothing, ic, bcp)

prob = HeatProblem(mesh, funs, paras)
sol = solve(prob, FDM())

@test sol.flow[1, 1, :, :] == sin.(2π * mesh[2, 1, :, :]) .^ 2
@test isapprox(maximum(sol.flow[1, end, :]), 0.0, atol = 1e-4)

# TODO find the reason why the result is oscillating, Runga's phenomenon ?

## test NSHeatProblem solved by NSFDM
mesh = zeros(3, nt, nx, ny)
for i in 1:nt, j in 1:nx, k in 1:ny
    mesh[1, i, j, k] = (i - 1) * lt / (nt - 1)
    mesh[2, i, j, k] = (j - 1) * lx / (nx - 1)
    mesh[3, i, j, k] = (k - 1) * (lye - lys) / (ny - 1) + lys
end

ic = zeros(2, nx, ny)
for i in 1:nx
    ic[1, i, :] .= sin(2π * xspan[i])^2
    ic[2, i, :] .= 4π * cos(2π * xspan[i]) * sin(2π * xspan[i])
end
fp = ones(1, nt, nx, ny)
bcp = zeros(2, nt, ny)
paras = (fp, nothing, ic, bcp)

prob = NSHeatProblem(mesh, Ny, funs, paras)
NSsol = solve(prob, NSFDM())

@test isapprox(NSsol.flow[1, 1, :, 1], sin.(2π .* xspan) .^ 2, atol = 1e-9)
@test isapprox(maximum(sol.flow[:, :, :] - NSsol.flow[:, :, :, 1]), 0.0, atol = 1e-9)
