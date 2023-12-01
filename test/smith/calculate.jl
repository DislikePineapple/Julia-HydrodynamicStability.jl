include("smith_equation.jl")

function eigenvalue!(fs::FreeStream, wa::Wave)
    prob = NonlinearProblem(dispersion, wa.α, [fs, wa])
    wa.α = solve(
        prob,
        Muller();
        abstol = 1e-9,
        maxiters = 1000,
        mode = "Sptial",
        # showiters = true,
    ).t
end

function eigenvalue_fix_point!(fs::FreeStream, wa::Wave)
    @unpack ω, β = wa

    ωfind = 0.5
    βfind = 0.5

    wa.α = 0.08 - 0.015im
    wa.ω = ωfind

    step = β >= βfind ? 0.01 : -0.01
    βs = βfind:step:β
    for β in βs
        wa.β = β
        eigenvalue!(fs, wa)
    end

    step = ω >= ωfind ? 0.01 : -0.01
    ωs = ωfind:step:ω
    for ω in ωs
        wa.ω = ω
        eigenvalue!(fs, wa)
    end
end

function eigenvalue_along_ω(fs::FreeStream, ω::AbstractArray, β::Number)
    ω[begin] < ω[end] && reverse!(ω)
    α = zeros(ComplexF64, length(ω))
    wa = Wave(ω[begin], Inf, β)
    for i in eachindex(ω)
        wa.ω = ω[i]
        if i == 1
            eigenvalue_fix_point!(fs, wa)
        elseif i == 2
            wa.α = α[i-1]
        else
            wa.α = 2 * α[i-1] - α[i-2]
        end
        α[i] = eigenvalue!(fs, wa)
    end
    return α
end
