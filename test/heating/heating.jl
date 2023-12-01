using HydrodynamicStability, Test, FFTW, LinearAlgebra

import Integrals
import SpecialFunctions: airyai, airyprime
import UnPack: @unpack

using Plots
plotlyjs()

include("type.jl")
include("carculate.jl")

fs = FreeStream(3, 273.15, Inf)
hs = Heating(0, 0.5)
wa = Wave(0.5, 0.08 - 0.013im, 0.5)

x = [1]
y = [Inf]
z = range(-π, π, 17)
z = z[1:16] / hs.β
gr = Grid(x, y, z)

eigenvalue!(fs, wa, hs, gr)
