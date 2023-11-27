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
Cprime_over_C(T, Tprime, Te) = 1 / 2 * Tprime / T * (Se / Te - T) / (T + Se / Te)
D(T, Te) = C(T, Te) * (1 / 2 / T - 1 / (T + Se / Te))
Dprime(T, Tprime, Te) =
    -Tprime * (Se + Te * T) * (Se^2 + 6 * Se * Te * T - 3 * Te^2 * T^2) /
    (4 * T^(3 / 2) * (Se + Te * T)^3)

function similarity!(du, u, p, t) # u[1] = f, u[2] = f', u[3] = g, u[4] = f'', u[5] = g'; f' = U, g = T , 
    @unpack Ma, Te = p

    du[1] = u[2]
    du[2] = u[4]
    du[3] = u[5]
    du[4] = -1 / 2 / C(u[3], Te) * u[1] * u[4] - Cprime_over_C(u[3], u[5], Te) * u[4]
    du[5] =
        -Pr / 2 / C(u[3], Te) * u[1] * u[5] - Cprime_over_C(u[3], u[5], Te) * u[5] -
        Pr * (Γ - 1) * Ma^2 * u[4]^2
end

function similarity_bc!(residual, u, p, t)
    Ma = p.Ma
    Tw = 1 + sqrt(Pr) * (Γ - 1) * Ma^2 / 2

    residual[1] = u[begin][1]       # F(0) = 0
    residual[2] = u[begin][2]       # U(0) = 0
    residual[3] = u[begin][3] - Tw  # T(0) = Tw
    residual[4] = u[end][2] - 1     # U → 1 as y → ∞
    residual[5] = u[end][3] - 1     # T → 1 as y → ∞
end
