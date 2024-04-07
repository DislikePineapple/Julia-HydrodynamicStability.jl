using HydrodynamicStability, NumericalIntegration, Plots, Revise
plotlyjs()
import UnPack: @unpack

include("type.jl")
include("equations.jl");
# parameter setting
Ma = 4.5
Te = 226.5
Re₁ = 5535304.0 # The unit Reynolds number

Re = 10456.106799816920  # The entrance Reynolds number
x₀ = 0.2 # The dimensional position for entrance Re number
δ₀ = Re / Re₁ # The dimensional displacement thickness at the entrance
δ = δ₀ / sqrt(x₀ / (Re₁ * x₀)) # The dimensional displacement thickness at dimensionless length L
fs = FreeStream(Ma, Te, Re₁)

grid, flow = compressibleBlasius(fs; δ = δ);

## 插值 生成新的网格和流动 --------------------------
using Interpolations

x = grid[1, :, 1]
ny_new = 401
grid_new = zeros(2, length(x), ny_new)
flow_new = zeros(4, length(x), ny_new)

δ_99 = findfirst(u -> u >= 0.99, flow[2, end, :])
y_mid = grid[2, end, δ_99]
y_max = grid[2, end, δ_99] * 10

a = y_mid * y_max / (y_max - 2.0 * y_mid)
b = 1 + 2 * a / y_max
y_fun(y, a, b) = a * (1 + y) / (b - y)
y_range = range(-1, 1, ny_new)

grid_new[1, :, :] = grid[1, :, 1:ny_new]
isapprox(grid_new[1, :, :], grid[1, :, 1:ny_new])
for i in eachindex(x)
    grid_new[2, i, :] = y_fun.(y_range, a, b)
end

for i in eachindex(x)
    j_max = findlast(grid_new[2, i, :] .< grid[2, i, end])
    for n = 1:4
        interp_linear = linear_interpolation(grid[2, i, :], flow[n, i, :])
        flow_new[n, i, 1:j_max] = interp_linear(grid_new[2, i, 1:j_max])
        for j = j_max+1:ny_new
            flow_new[n, i, j] = flow_new[n, i, j_max]
        end
    end
end

# ∂(ρu)∂x + ∂(ρv)∂y = 0
ρuₓ, v = [zeros(size(flow_new[2, :, :])) for _ = 1:2]
# ∂(ρu)∂x
for j in eachindex(grid_new[2, 1, :])
    ρuₓ[:, j] = HydrodynamicStability.central_difference(
        flow_new[1, :, j] .* flow_new[2, :, j],
        grid_new[1, :, j],
    )
end
# ∂(ρv)∂y
for i in eachindex(grid[1, :, 1])
    v[i, :] =
        HydrodynamicStability.simpsons_integral(-ρuₓ[i, :], grid_new[2, i, :]) ./
        flow_new[1, i, :]
end

# plot for check v ---------------------------------
plot(grid_new[1, :, 1], v[:, end])
plot!(grid_new[1, :, 1], flow[3, :, end])

##* write as Plot3D files --------------------------
using FortranFiles

file = FortranFile("data/grid.xyz", "w")
write(file, Int32(1))
write(file, Int32(length(x)), Int32(ny_new))
write(file, grid_new[1, :, :], grid_new[2, :, :])
close(file)

file = FortranFile("data/flow.fun", "w")
write(file, Int32(1))
write(file, Int32(length(x)), Int32(ny_new), Int32(4))
write(file, flow_new[1, :, :], flow_new[2, :, :], flow_new[3, :, :], flow_new[4, :, :])
close(file)
