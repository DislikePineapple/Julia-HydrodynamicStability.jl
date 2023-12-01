const λ₀ = 0.3321
const Pr = 0.72
const Γ = 1.4
const Se = 110.4

function baseflow(fs::FreeStream)
    @unpack Ma, Te = fs
    t_B = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2
    ρ_B = t_B^(-1)
    C_B = t_B^(1 / 2) * (1 + Se / Te) / (t_B + Se / Te)

    return t_B, ρ_B, C_B
end

function baseflow(fs::FreeStream, hs::Heating, gr::Grid)
    @unpack Ma, Te = fs
    @unpack h, f, s = hs
    @unpack x, z = gr

    t_B0 = 1 + 0.5 * Pr^(1 / 2) * (Γ - 1) * Ma^2

    t_w = t_B0 .+ h * t_B0 .* f.(x) .* s.(z)
    t_wz = central_difference(t_w, z)

    ρ_w = t_w .^ (-1)
    ρ_wz = central_difference(ρ_w, z)

    λ_u = λ₀ * t_B0 ./ t_w
    λ_uz = central_difference(λ_u, z)

    μ_w = t_w .^ (3 / 2) .* (1 + Se / Te) ./ (t_w .+ Se / Te)
    μ_wz = central_difference(μ_w, z)

    d_w = 1.5 ./ t_w - 1.0 ./ (t_w .+ Se / Te)

    return t_w, t_wz, λ_u, λ_uz, μ_w, μ_wz, ρ_w, ρ_wz, d_w
end

C1 = [
    0 0 -1.5 2 -0.5
    0 -1/3 -0.5 1 -1/6
    1/12 -2/3 0 2/3 -1/12
    1/6 -1 0.5 1/3 0
    1/2 -2 1.5 0 0
]

function central_difference(f, x)
    df = zeros(eltype(f), length(f))
    df[1] = (C1[1, 3] * f[1] + C1[1, 4] * f[2] + C1[1, 5] * f[3]) / (x[2] - x[1])
    df[2] =
        (C1[2, 2] * f[1] + C1[2, 3] * f[2] + C1[2, 4] * f[3] + C1[2, 5] * f[4]) /
        (x[3] - x[1]) * 2
    for j = 3:length(f)-2
        df[j] =
            (
                C1[3, 1] * f[j-2] +
                C1[3, 2] * f[j-1] +
                C1[3, 3] * f[j] +
                C1[3, 4] * f[j+1] +
                C1[3, 5] * f[j+2]
            ) / (x[j+1] - x[j-1]) * 2
    end
    df[end-1] =
        (
            C1[4, 1] * f[end-3] +
            C1[4, 2] * f[end-2] +
            C1[4, 3] * f[end-1] +
            C1[4, 4] * f[end]
        ) / (x[end] - x[end-2]) * 2
    df[end] =
        (C1[5, 1] * f[end-2] + C1[5, 2] * f[end-1] + C1[5, 3] * f[end]) /
        (x[end] - x[end-1])
    return df
end
