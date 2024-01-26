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
function C(T, Te)
    T^(1 / 2) * (1 + Se / Te) / (T + Se / Te)
end
Cprime_over_C(T, Tprime, Te) = Tprime * (Se / Te - T) / (T + Se / Te) / T / 2

"""
    similarity!(du, u, p, t)

The compressible blasius equations:

```math
    -\frac{1}{2}ff''=(Cf'')',
    -\frac{1}{2}fT'=\frac{1}{P_r}(CT')'+(γ-1)Ma^2C(f'')^2.
```

"""
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

"""
    similarity_bc!(residual, u, p, t)

The boundary conditions for compressible blasius equations, which means:
    
```math
    f(0) = 0, f'(0) = 0, f' → 1 as y → ∞,
    ∫_0^∞ T(y)dy = 0, T(0) = T_w, T → 1 as y → ∞.
```
    
"""
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
