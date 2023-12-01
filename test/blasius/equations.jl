"""
## Incompressible blasius equations
``` 
f'''+1/2ff''=0
```
"""
function blasius!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
end

function blasius_bc!(residual, u, p, t)
    residual[1] = u[begin][1]    # f(0)  = 0
    residual[2] = u[begin][2]    # f'(0) = 0
    residual[3] = u[end][2] - 1  # f' → 1 as y → ∞
end

C(T, Te) = T^(1 / 2) * (1 + Se / Te) / (T + Se / Te)
Cprime_over_C(T, Tprime, Te) = Tprime * (Se / Te - T) / (T + Se / Te) / T / 2
# Cprime_over_C(T, Tprime, Te) = Tprime * (0.5 / T - 1 / (T + Se / Te))

# C(T, Te) = 1.5
# Cprime_over_C(T, Tprime, Te) = 0

function similarity!(du, u, p, t) # u[1] = f, u[2] = f', u[3] = g, u[4] = f'', u[6] = g'; f' = U, g = T , 
    @unpack Ma, Te = p

    du[1] = u[2]
    du[2] = u[3]
    du[3] = -u[1] * u[3] / C(u[5], Te) / 2 - Cprime_over_C(u[5], u[6], Te) * u[3]
    du[4] = u[5]
    du[5] = u[6]
    du[6] =
        -Pr * u[1] * u[6] / C(u[5], Te) / 2 - Cprime_over_C(u[5], u[6], Te) * u[6] -
        Pr * (Γ - 1) * Ma^2 * u[3]^2
end

function similarity_bc!(residual, u, p, t)
    Ma = p.Ma
    Tw = 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2

    residual[1] = u[begin][1]       # F(0) = 0
    residual[2] = u[begin][2]       # U(0) = 0
    residual[3] = u[end][2] - 1     # U → 1 as y → ∞
    residual[4] = u[begin][4]       # ∫T(0) = 0
    residual[5] = u[begin][5] - Tw  # T(0) = Tw
    residual[6] = u[end][5] - 1     # T → 1 as y → ∞
end
