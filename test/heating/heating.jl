using HydrodynamicStability, FFTW, LinearAlgebra, ProgressMeter

import Integrals
import SpecialFunctions: airyai, airyaiprime
import UnPack: @unpack

include("types.jl")

include("equations/baseflow.jl")
include("equations/function.jl")

include("equations/smith_equation.jl")
include("equations/heating_equation.jl")

function eigenvalue!(fs::FreeStream, wa::Wave; kwarg...)
    prob = NonlinearProblem(smith, wa.α, [fs, wa])
    wa.α = solve(
        prob,
        Muller();
        kwarg...,
        # abstol = 1e-9,
        # maxiters = 1000,
        # mode = "Sptial",
        # showiters = true,
    ).t
end

function eigenvalue!(fs::FreeStream, wa::Wave, hs::Heating; kwarg...)
    prob = NonlinearProblem(eigenproblem, wa.α, [fs, wa, hs])
    wa.α = solve(
        prob,
        Muller();
        # abstol = 1e-9,
        # maxiters = 1000,
        # showiters = true,
        # δ = 1e-8im,
        Return = "eigenvalue",
        kwarg...,
    ).t
end

## two-dimentional functions for carculation --------------------------------------------------------------------

function eigenvalue_fix_point!(fs::FreeStream, wa::Wave; kwarg...)
    @unpack ω, β = wa

    ωfind = 0.5
    βfind = 0.5

    wa.α = 0.08 - 0.015im
    wa.ω = ωfind

    step = β >= βfind ? 0.01 : -0.01
    βs = βfind:step:β
    for β in βs
        wa.β = β
        eigenvalue!(fs, wa; kwarg...)
    end

    step = ω >= ωfind ? 0.01 : -0.01
    ωs = ωfind:step:ω
    for ω in ωs
        wa.ω = ω
        eigenvalue!(fs, wa; kwarg...)
    end
    return wa.α
end

function eigenvalue_along_ω(
    fs::FreeStream,
    wa::Wave,
    ωspan::AbstractArray;
    initial = true,
    kwarg...,
)
    ωspan[begin] < ωspan[end] && reverse!(ωspan)
    αspan = zeros(ComplexF64, length(ωspan))
    for i in eachindex(ωspan)
        wa.ω = ωspan[i]
        if i == 1
            initial && eigenvalue_fix_point!(fs, wa; kwarg...)
        elseif i == 2
            wa.α = αspan[i-1]
        else
            wa.α = 2 * αspan[i-1] - αspan[i-2]
        end
        αspan[i] = eigenvalue!(fs, wa; kwarg...)
    end
    # println("test")
    return αspan
end


## 3D --------------------------------------------------------------------

function eigenvalue_fix_point!(fs::FreeStream, wa::Wave, hs::Heating; kwarg...)
    xspan = 0:0.01:hs.x
    eigenvalue_along_x(fs, wa, hs, xspan; kwarg...)[end]
end

function eigenvalue_along_x(
    fs::FreeStream,
    wa::Wave,
    hs::Heating,
    xspan::AbstractArray;
    initial = true,
    kwarg...,
)
    initial && (wa.α = eigenvalue_fix_point!(fs, wa; kwarg...))
    αspan = zeros(ComplexF64, length(xspan))
    for i in eachindex(xspan)
        if i == 1
        elseif i == 2
            wa.α = αspan[i-1]
        else
            wa.α = 2 * αspan[i-1] - αspan[i-2]
        end
        hs.x = xspan[i]
        αspan[i] = eigenvalue!(fs, wa, hs; kwarg...)
    end
    return αspan
end

function eigenvalue_along_x(
    fs::FreeStream,
    wa::Wave,
    hs::Heating,
    xspan::AbstractArray,
    ωspan::AbstractArray;
    kwarg...,
)
    αplane = zeros(ComplexF64, length(ωspan), length(xspan))
    αplane[:, 1] = eigenvalue_along_ω(fs, wa, ωspan; kwarg...)

    @showprogress 1 "Eigenvalue for x-omega plane" for i in eachindex(ωspan)
        wa.ω = ωspan[i]
        wa.α = αplane[i, 1]
        αplane[i, :] = eigenvalue_along_x(fs, wa, hs, xspan; initial = false, kwarg...)
    end
    return αplane
end

function eigenvalue_along_ω(
    fs::FreeStream,
    wa::Wave,
    hs::Heating,
    ωspan::AbstractArray;
    initial = true,
    kwarg...,
)
    ωspan[begin] < ωspan[end] && reverse!(ωspan)
    αspan = zeros(ComplexF64, length(ωspan))
    @showprogress 1 "Eigenvalue for omega span" for i in eachindex(ωspan)
        wa.ω = ωspan[i]
        if i == 1
            initial && eigenvalue_fix_point!(fs, wa, hs; kwarg...)
        elseif i == 2
            wa.α = αspan[i-1]
        else
            wa.α = 2 * αspan[i-1] - αspan[i-2]
        end
        αspan[i] = eigenvalue!(fs, wa, hs; kwarg...)
    end
    return αspan
end

function eigenvalue_along_ω(
    fs::FreeStream,
    wa::Wave,
    hs::Heating,
    ωspan::AbstractArray,
    xspan::AbstractArray;
    initial = true,
    kwarg...,
)
    ωspan[begin] < ωspan[end] && reverse!(ωspan)
    αplane = zeros(ComplexF64, length(ωspan), length(xspan))
    αplane[:, 1] = eigenvalue_along_ω(fs, wa, ωspan; kwarg...)
    @showprogress 1 "Eigenvalue for x-omega plane" for j in eachindex(xspan)
        hs.x = xspan[j]
        if j == 1
            initial && eigenvalue_fix_point!(fs, wa; kwarg...)
            continue
        elseif j == 2
            wa.α = αplane[1, j-1]
        else
            wa.α = 2 * αplane[1, j-1] - αplane[1, j-2]
        end
        αplane[:, j] = eigenvalue_along_ω(fs, wa, hs, ωspan; kwarg...)
    end
    return αplane
end
