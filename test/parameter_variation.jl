using HydrodynamicStability, Test, SpecialFunctions

"""
    L _Fun!(D, F, p, t)

```math
Ai'(x) = Bi(x) ∫_{x}^{∞} Ai(y) Ai(y) dx + Ai(x) ∫_{0}^{x} Ai(y) Bi(y) dx
```

"""

yspan = range(-10, 10, 1001)
fun = airyai.(yspan)

sol = parameter_variation(fun, collect(yspan))
sol = -sol[findfirst(x -> x == 0, yspan)] / airyai(0) * airyai.(yspan) + sol
solm = -airyaiprime(0) / airyai(0) * airyai.(yspan) + airyaiprime.(yspan)

@test isapprox(sol, solm, atol = 1e-5)

##* Nonlinear Equation ----------------------------------------------

function functionR(Pr, ζspan)
    function L_Fun!(A, D, F, t, p)
        A[1, 1] = 1
        A[2, 2] = 1

        D[1, 2] = -1
        D[2, 1] = -t

        F[1] = 0
        F[2] = 1
    end
    function L_bc!(M, u)
        M[2, :] .= 0
        M[2, 1] = 1
        M[end, :] .= 0
        M[end, end - 1] = 1

        u[2] = 0
        u[end] = 0
    end

    u₀ = complex([0.0, 0.0])
    prob = BVProblem(L_Fun!, L_bc!, u₀, ζspan)
    ret = solve(prob, FDM(); ny = length(ζspan))
    sol = flatten_vector(ret.u)

    L = sol[1, :]

    ζspan_Pr = Pr^(1 / 3) .* ζspan

    sol_LT = parameter_variation(L, collect(ζspan_Pr))
    sol_LT = -sol_LT[findfirst(x -> x == ζspan_Pr[1], ζspan_Pr)] / airyai(ζspan_Pr[1]) .*
             airyai.(ζspan_Pr) + sol_LT

    solLTDoublePrime = central_difference(
        central_difference(sol_LT, ζspan_Pr; accurency = 4), ζspan_Pr; accurency = 4)
    sol_RT = parameter_variation(solLTDoublePrime, collect(ζspan))
    D_sol_RT = central_difference(sol_RT, ζspan; accurency = 4)
    sol_RT = -D_sol_RT[findfirst(x -> x == ζspan[1], ζspan)] / airyaiprime(ζspan[1]) .*
             airyai.(ζspan) + sol_RT
    Int_RT = simpsons_integral(sol_RT, ζspan; accurency = 3)

    return sol, sol_LT, Int_RT
end

yspan = range(1, 20 + 1, 101) * exp(π / 6 * 1im)
y0 = exp(π / 6 * 1im)
y0N = findfirst(x -> x == y0, yspan)

sol_L, sol_LT, Int_RT = functionR(1, yspan)
sol_L_Pr, sol_LT_Pr, Int_RT_Pr = functionR(0.72, yspan)

L = sol_L[1, :]
Lprime = sol_L[2, :]

sol_LT_A = -Lprime[y0N] / airyai(y0) .* airyai.(yspan) + Lprime

κ = simpsons_integral(airyai.(yspan[y0N:end]),
    yspan[y0N:end]; accurency = 3)
Phi00 = Lprime[y0N] / 3 / airyai(y0) *
        (airyaiprime(y0)) - 3 / 4 + Lprime[y0N] / 12 * y0^2
Int_RT_A = Phi00 / airyaiprime(y0) .* κ +
           Lprime[y0N] / 3 / airyai(y0) .* κ -
           Lprime[y0N] / 3 / airyai(y0) .*
           yspan[y0N:end] .*
           airyai.(yspan[y0N:end]) +
           L[y0N:end] ./ 4 +
           yspan[y0N:end] .*
           Lprime[y0N:end] ./ 4 .+
           y0 * Lprime[y0N] / 12

@test isapprox(sol_LT_A[1:51], sol_LT[1:51], atol = 1e-2)
@test isapprox(sum(abs2, Int_RT_A[1:51] - Int_RT[1:51]), 0.0, atol = 1e-2)
