using HydrodynamicStability, Test

import Integrals
import SpecialFunctions: airyai, airyaiprime
import UnPack: @unpack

using Plots
plotlyjs()

# include("type.jl")
# include("calculate.jl")

fs = FreeStream(3, 273.15, Inf)
ω = collect(0.0:0.001:2)
β = 0.5
# α = eigenvalue_along_ω(fs, ω, β)

wa = Wave(0.5, 0.08 - 0.015im, 0.5)
eigenvalue!(fs, wa)
# 0.07905835533381554 - 0.013615680374057752im

# plt = plot(
#     1,
#     xaxis = ("𝜔", (0, 1), font(15, "Times")),
#     yaxis = ("-𝛼ᵢ", (-0.01, 0.03), font(15, "Times")),
#     framestyle = :box,
#     fontfamily = "Times",
#     legendfontsize = 12,
# )

# plot(plt, ω, -imag(α), label = nothing, w = 4, lc = :black)
# plt_imag = plot!(
#     xaxis = ("𝜔", (0, 1), font(15, "Times")),
#     yaxis = ("-𝛼ᵢ", (-0.01, 0.03), font(15, "Times")),
#     framestyle = :box,
#     fontfamily = "Times",
#     legendfontsize = 12,
# )
# wa = Wave(0.5, 0.08 - 0.01im, 0.4)
# eigenvalue_fix_point!(fs, wa)
# wa.α
# eigenvalue!(fs, wa)
