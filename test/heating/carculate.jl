include("heating_equation.jl")

function eigenvalue!(fs::FreeStream, wa::Wave, hs::Heating, gr::Grid)
    prob = NonlinearProblem(eigenproblem, wa.α, [fs, wa, hs, gr])
    wa.α =
        solve(
            prob,
            Muller();
            abstol = 1e-9,
            maxiters = 1000,
            showiters = true,
            δ = 0.0001im,
            Return = "eigenvalue",
        ).t
end
