using HydrodynamicStability, Test

import Integrals
import SpecialFunctions: airyai, airyaiprime
import UnPack: @unpack

include("type.jl")
include("calculate.jl")

fs = FreeStream(3, 250, Inf)
ω = collect(0.01:0.01:10)
β = 2
α = eigenvalue_along_ω(fs, ω, β)

using Plots
plotlyjs()
plot(ω, -imag(α), label = nothing, w = 4, lc = :black)
plt_imag = plot!(
    xaxis = ("𝜔", (0, 10), font(15, "Times")),
    yaxis = ("-𝛼ᵢ", (-0.01, 0.05), font(15, "Times")),
    framestyle = :box,
    fontfamily = "Times",
    legendfontsize = 12,
)
