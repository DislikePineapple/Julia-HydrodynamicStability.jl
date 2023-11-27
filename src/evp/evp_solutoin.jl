struct EVPSolution{V,E,P,A,RT} <: AbstractEVPSolution
    ev::V
    ef::E

    prob::P
    alg::A

    retcode::RT
end

EVPSolution(ev) = EVPSolution(ev, nothing)
EVPSolution(ev, ef) = EVPSolution(ev, ef, nothing, nothing)
EVPSolution(ev, ef, prob, alg) = EVPSolution(ev, ef, prob, alg, nothing)
