using HydrodynamicStability, NumericalIntegration, Plots
plotlyjs()
import UnPack: @unpack

include("type.jl")
include("equations.jl")
;
## parameter setting
Ma = 3
Te = 216.7
Re1 = 5535304.0
Re = 10456.106799816920
fs = FreeStream(Ma, Te, Re)

u₀ = [0, 0, 0.3779, 0, 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2, 0]
tspan = (0.0, 100.0)

# solve the compressible blasius equations in the coodinate η = 1 / sqrt(x * Re) ∫_0^y ρ(η) dη
bvp = BVProblem(similarity!, similarity_bc!, u₀, tspan, fs)
sol = solve(bvp, Shooting(), dy = 0.05);

## postporcessing
# define the domain and give the streamwise grid
x = 0.75:0.005:4
grid = zeros(2, length(x), length(sol.y)) # grid[1,:,:] = x, grid[2,:,:] = y
for j = 1:length(sol.y)
    grid[1, :, j] = x
end

# define the flow field ρ, u, v, T
flow = zeros(4, length(x), length(sol.y))
y, F = [zeros(length(sol.y)) for i = 1:2]
for i in eachindex(x), j in eachindex(sol.y)
    flow[1, i, j] = 1 / sol.u[j][5]
    flow[2, i, j] = sol.u[j][2]
    flow[4, i, j] = sol.u[j][5]
    F[j] = sol.u[j][1]
    y[j] = sol.u[j][4]
end

# carculate the entrance momentum thickness
ρu = flow[1, 1, :] .* flow[3, 1, :]
ratio = Re / Re1 / integrate(y, 1 .- ρu) / sqrt(0.02)
# carculate normalized grid
for i in eachindex(x)
    grid[2, i, :] = sqrt(x[i]) * y * ratio
end

# solve normal velocity v
for i in eachindex(x)
    flow[3, i, :] = -1 / sqrt(x[i] / 0.02 * Re) * (F .* flow[4, 1, :] - y .* flow[2, 1, :])
end

## 插值 生成新的网格和流动
using Interpolations

ny_new = 400
y_j = zeros(ny_new);
grid_new = zeros(3, length(x), ny_new)

grid_new[1, :, 1:ny_new] = grid[1, :, 1:ny_new]
flow_new = zeros(4, length(x), ny_new)
for j = 0:ny_new-1
    y_j[ny_new-j] = 1 - (2.0 / (ny_new - 1)) * j   #cos(j*pi/(ny_new-1))
end;
y_new(ŷ, a, b) = a * (1 + ŷ) / (b - ŷ);
y_i = grid[2, end, 120];
y_max = grid[2, end, end];

for i = 1:length(x)
    a = y_i * y_max / (y_max - 2.0 * y_i)
    b = 1 + 2 * a / y_max
    grid_new[2, i, :] = y_new.(y_j, a, b)

    ỹ = grid_new[2, i, :]
    interp_linear = LinearInterpolation(grid[2, i, :], flow[1, i, :])
    j_max  = length(ỹ[isless.(ỹ, grid[2, i, end])])
    flow_new[1, i, 1:j_max] = interp_linear(grid_new[2, i, 1:j_max])
    interp_linear = LinearInterpolation(grid[2, i, :], flow[2, i, :])
    flow_new[2, i, 1:j_max] = interp_linear(grid_new[2, i, 1:j_max])
    interp_linear = LinearInterpolation(grid[2, i, :], flow[3, i, :])
    flow_new[3, i, 1:j_max] = interp_linear(grid_new[2, i, 1:j_max])
    interp_linear = LinearInterpolation(grid[2, i, :], flow[4, i, :])
    flow_new[4, i, 1:j_max] = interp_linear(grid_new[2, i, 1:j_max])
    j_max = j_max + 1
    for j = j_max:ny_new
        flow_new[1, i, j] = flow_new[1, i, j_max-1]
        flow_new[2, i, j] = flow_new[2, i, j_max-1]
        flow_new[3, i, j] = flow_new[3, i, j_max-1]
        flow_new[4, i, j] = flow_new[4, i, j_max-1]
    end
end

# plot for check
plot(flow_new[3, 100, 1:150], grid_new[2, 1, 1:150])

## write as Plot3D files
using FortranFiles

f = FortranFile("data/grid.xyz", "w")
write(f, Int32(1))
write(f, Int32(length(x)), Int32(ny_new))
write(f, grid_new[1, :, :], grid_new[2, :, :])
close(f)

f = FortranFile("data/flow.fun", "w")
write(f, Int32(1))
write(f, Int32(length(x)), Int32(ny_new), Int32(4))
write(f, flow_new[1, :, :], flow_new[2, :, :], flow_new[3, :, :], flow_new[4, :, :])
close(f)
