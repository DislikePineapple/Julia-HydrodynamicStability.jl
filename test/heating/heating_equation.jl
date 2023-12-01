include("baseflow.jl")

function eigenproblem(α::ComplexF64, p; Return)
    fs = p[1]
    wa = p[2]
    hs = p[3]
    gr = p[4]

    β = hs.β
    q = wa.β / β

    wa.α = α

    Nz = length(gr.z)
    NM = Int(Nz / 2 - 1)
    M = zeros(ComplexF64, NM, NM)

    # function eigenproblem(α)
    L, Q = functionLQ(fs, wa, hs, gr)

    for j = 1:NM
        nnum = j - Nz / 4 - 1
        γs = γ(fs.Ma, α, β, q, nnum)

        for n = 1:NM
            M[n, j] =
                im * β * (q + nnum) * L[Int(Nz / 2 + 1 + n - j)] +
                1 / α^2 * γs * Q[Int(Nz / 2 + 1 + n - j)]
        end
        M[j, j] += α^2 + β^2 * (q + nnum)^2
    end

    Return == "eigen_matrix" && return M
    Return == "eigenvalue" && return det(M)
    # M = lu(M1 + M2)
    # Return == "dispersion" && return M.U[end, end]
    error("return message incrroct in eigenproblem! choose eigen_matrix or dispersion.")
    # end
    # return eigenproblem
end

airyai(t, p) = airyai(t)
κ(ζ₀, ζ) = Integrals.solve(Integrals.IntegralProblem(airyai, ζ₀, ζ), Integrals.QuadGKJL()).u

function γ(Ma, α, β, q, n)
    δ = β^2 * (q + n)^2 + α^2 * (1 - Ma^2)
    if real(δ) > 0
        γs = √δ
    else
        γs = im * √ - δ
    end
    return γs
end

function functionLQ(fs, wa, hs, gr)
    @unpack ω, α = wa
    Nz = length(z)

    _, t_wz, λ_u, λ_uz, μ_w, μ_wz, ρ_w, ρ_wz, d_w = baseflow(fs, hs, gr)

    ζ₀ = -(im * α .* λ_u .* ρ_w ./ μ_w) .^ (1 / 3) .* ω ./ (α .* λ_u)

    Ai₀ = airyai.(ζ₀)
    Aiprime₀ = airyprime.(ζ₀)

    κ₀ = 1 / 3 .- κ.(zero(eltype(ζ₀)), ζ₀)

    L1 =
        (0.25 * μ_wz ./ μ_w .- 0.25 * ρ_wz ./ ρ_w .+ 0.5 * λ_uz ./ λ_u) ./ Ai₀ .*
        (ζ₀ .^ 2 .* κ₀ .+ ζ₀ .* Aiprime₀)
    L2 = (0.25 * μ_wz ./ μ_w .+ 0.75 * ρ_wz ./ ρ_w .+ 1.5 * λ_uz ./ λ_u)
    L3 =
        (
            1 / 12 ./ Ai₀ .* (ζ₀ .^ 2 .* κ₀ .+ ζ₀ .* Aiprime₀) .+ 0.75 .+
            2 / 3 * κ₀ .* Aiprime₀ ./ Ai₀ .^ 2
        ) .* (ρ_w .- d_w) .* t_wz
    L4 = (1 .+ Aiprime₀ .* κ₀ ./ Ai₀ .^ 2) .* d_w .* t_wz
    L = L1 .+ L2 .+ L3 .+ L4
    Q = (im .* α .* λ_u) .^ (5 / 3) .* ρ_w .^ (2 / 3) .* μ_w .^ (1 / 3) .* Aiprime₀ ./ κ₀

    L = fftshift(fft(L) ./ Nz)
    Q = fftshift(fft(Q) ./ Nz)
    return L, Q
end
