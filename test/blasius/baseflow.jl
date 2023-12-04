function baseflow(fs::FreeStream)
    @unpack Ma, Te = fs
    t_B = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2
    ρ_B = t_B^(-1)
    C_B = t_B^(1 / 2) * (1 + Se / Te) / (t_B + Se / Te)

    return t_B, ρ_B, C_B
end
