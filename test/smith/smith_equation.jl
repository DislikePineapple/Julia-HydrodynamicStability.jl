# include("baseflow.jl")
const λ₀ = 0.3321
const Pr = 0.72
const Γ = 1.4
const Se = 110.4

airyai(t, p) = airyai(t)
κ(ζ₀, ζ) = Integrals.solve(Integrals.IntegralProblem(airyai, ζ₀, ζ), Integrals.QuadGKJL()).u

function dispersion(val::ComplexF64, p; mode)
    Ma = p[1].Ma
    @unpack ω, α, β = p[2]
    _, ρ_B, μ_B = baseflow(p[1])

    if mode == "Temporal"
        ω = val
    elseif mode == "Sptial"
        α = val
    else
        error("The mode is not supported in dispersion! Choose Temporal or Sptial.")
    end

    γ = (β^2 - α^2 * (Ma^2 - 1))^(1 / 2)
    ζ₀ = -(im * α * λ₀ * ρ_B / μ_B)^(1 / 3) * ω / (α * λ₀)
    residual =
        α^2 * (α^2 + β^2) +
        (im * α * λ₀)^(5 / 3) * ρ_B^(2 / 3) * μ_B^(1 / 3) * γ * airyaiprime(ζ₀) /
        (1 / 3 - κ(zero(eltype(ζ₀)), ζ₀))

    # γ = (β^2 / (Ma^2 - 1) - α^2)^(1 / 2)
    # ζ₀ = -(im * α)^(1 / 3) * ω / α
    # residual =
    #     (im * α)^(1 / 3) * (α^2 + β^2) -
    #     γ * airyaiprime(ζ₀) / (1 / 3 - κ(zero(eltype(ζ₀)), ζ₀))

end

function baseflow(fs::FreeStream)
    @unpack Ma, Te = fs
    t_B = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2
    ρ_B = t_B^(-1)
    μ_B = t_B^(3 / 2) * (1 + Se / Te) / (t_B + Se / Te)

    return t_B, ρ_B, μ_B
end
