"""
    blasius!(du, u, p, t)

The incompressible blasius equations:

```math
    f'''+\frac{1}{2}ff''=0, 
```

which can be written as a first order system:

```math
    (f)'=f',
    (f')'=f'',
    (f'')'=-\frac{1}{2}ff''.
```

"""
function blasius!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
end

"""
    blasiusBC!(residual, u, p, t)

The boundary conditions for blasius equations, which means:

```math
    f(0) = 0, f'(0) = 0, f' → 1 as y → ∞.
```

addtion
"""
function blasius_bc!(residual, u, p, t)
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[end][2] - 1
end

"""
    C(T,Te)

The coefficient for Sutherland's viscosity law:

```math
    μ=CT=T^{3/2}\frac{1+S/T_e}{T+S/T_e}, 
```

"""
# C(T, Te) = T^(1 / 2) * (1 + Se / Te) / (T + Se / Te)
# Cprime_over_C(T, Tprime, Te) = Tprime * (Se / Te - T) / (T + Se / Te) / T / 2

# Chapman's law
C(T, Te) = 1
Cprime_over_C(T, Tprime, Te) = 0

"""
    similarity!(du, u, p, t)

The compressible blasius equations:

```math
    -\frac{1}{2}ff''=(Cf'')',
    -\frac{1}{2}fT'=\frac{1}{P_r}(CT')'+(γ-1)Ma^2C(f'')^2.
```

"""
function similarity!(du, u, fs, t) # u[1] = f, u[2] = f', u[3] = g, u[4] = f'', u[6] = g'; f' = U, g = T , 
    @unpack Ma, Te = fs

    du[1] = u[2]
    du[2] = u[3]
    du[3] = -u[1] * u[3] / C(u[5], Te) - Cprime_over_C(u[5], u[6], Te) * u[3]
    du[4] = u[5]
    du[5] = u[6]
    du[6] =
        -Pr * u[1] * u[6] / C(u[5], Te) - Cprime_over_C(u[5], u[6], Te) * u[6] -
        Pr * (γ - 1) * Ma^2 * u[3]^2
end

"""
    similarity_bc!(residual, u, p, t)

The boundary conditions for compressible blasius equations, which means:
    
```math
    f(0) = 0, f'(0) = 0, f' → 1 as y → ∞,
    ∫_0^∞ T(y)dy = 0, T(0) = T_w, T → 1 as y → ∞.
```
    
"""
function similarity_bc!(residual, u, fs, t)
    Tw = 1 + sqrt(Pr) * (γ - 1) * fs.Ma^2 / 2

    residual[1] = u[begin][1]       # F(0) = 0
    residual[2] = u[begin][2]       # U(0) = 0
    residual[3] = u[end][2] - 1     # U → 1 as y → ∞
    residual[4] = u[begin][4]       # ∫T(0) = 0
    residual[5] = u[begin][5] - Tw  # T(0) = Tw
    residual[6] = u[end][5] - 1     # T → 1 as y → ∞
end

function compressibleBlasius(fs::FreeStream; δ = nothing)
    u₀ = [0, 0, 0.3779, 0, 1 + sqrt(Pr) * (γ - 1) * fs.Ma^2 / 2, 0]
    # initial guess for the shooting method
    ηspan = (0.0, 100.0)

    # solve the compressible blasius equations in the coodinate η = 1 / sqrt(x * Re) ∫_0^y ρ(η) dη
    bvp = BVProblem(similarity!, similarity_bc!, u₀, ηspan, fs)
    sol = solve(bvp, Shooting(), dy = 0.05)

    ## postporcessing
    # define the domain and give the streamwise grid
    x = 0.75:0.005:4
    grid = zeros(2, length(x), length(sol.y))
    for j = 1:length(sol.y)
        grid[1, :, j] = x
    end

    # define the flow field ρ, u, v, T
    flow = zeros(4, length(x), length(sol.y))
    y, f = [zeros(length(sol.y)) for _ = 1:2]
    for j in eachindex(sol.y)
        flow[1, :, j] .= 1 / sol.u[j][5] # ρ
        flow[2, :, j] .= sol.u[j][2] # u
        flow[4, :, j] .= sol.u[j][5] # T
        f[j] = sol.u[j][1] # f
        y[j] = sol.u[j][4] # y = ∫_0^η ρ(η*) dη*
    end

    # carculate the entrance momentum thickness
    ρu = flow[1, 1, :] .* flow[2, 1, :]
    δ₉₉ = sol.y[findfirst(u -> u >= 0.99, flow[2, 1, :])]
    δd = integrate(sol.y, 1 .- ρu)
    δm = integrate(sol.y, ρu .* (1 .- ρu))

    if isnothing(δ) || isinf(fs.Re)
        # carculate normalized grid
        for i in eachindex(x)
            grid[2, i, :] = sqrt(2 * x[i]) * y
        end
        # solve normal velocity v
        for i in eachindex(x)
            flow[3, i, :] = -1 / sqrt(2 * x[i]) * (f .* flow[4, i, :] - y .* flow[2, i, :])
        end
    else
        ratio = δd / δ
        L = ratio^2
        @show ratio, L
        @show δ₉₉, δd, δm
        for i in eachindex(x)
            grid[2, i, :] = sqrt(2 * x[i] / fs.Re) * y * ratio
        end
        # solve normal velocity v
        for i in eachindex(x)
            flow[3, i, :] =
                -1 / sqrt(2 * x[i] / L * fs.Re) * (f .* flow[4, i, :] - y .* flow[2, i, :])
        end
    end
    return grid, flow
end
