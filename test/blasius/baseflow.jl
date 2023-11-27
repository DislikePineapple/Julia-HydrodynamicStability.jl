function baseflow(p)
    @unpack Ma, Te = p
    t_B = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2
    ρ_B = t_B^(-1)
    μ_B = t_B^(3 / 2) * (1 + Se / Te) / (t_B + Se / Te)

    return t_B, ρ_B, μ_B
end
