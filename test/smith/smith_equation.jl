# include("baseflow.jl")

airyai(t, p) = airyai(t)
κ(ζ₀, ζ) = Integrals.solve(Integrals.IntegralProblem(airyai, ζ₀, ζ), Integrals.QuadGKJL()).u

function dispersion(val, p; mode)
    Ma = p[1].Ma
    @unpack ω, α, β = p[2]

    if mode == "Temporal"
        ω = val
    elseif mode == "Sptial"
        α = val
    else
        error("The mode is not supported in dispersion! Choose Temporal or Sptial.")
    end

    γ = (β^2 / (Ma^2 - 1) - α^2)^(1 / 2)
    ζ₀ = -(im * α)^(1 / 3) * ω / α
    residual =
        (im * α)^(1 / 3) * (α^2 + β^2) -
        γ * airyaiprime(ζ₀) / (1 / 3 - κ(zero(eltype(ζ₀)), ζ₀))
end

function smith_ode(du, u, t, p)
    @unpack α, ω = p
    du[1] = u[2]
    du[2] = im * (α * t - ω) * u[1]
end
