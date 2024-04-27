struct NonlinearSolution{tType, P, A, tType2, RT} <: AbstractNonlinearSolution
    t::tType

    prob::P
    alg::A

    t_analytic::tType2
    retcode::RT
end

NonlinearSolution(t) = NonlinearSolution(t, nothing, nothing)
NonlinearSolution(t, prob, alg) = NonlinearSolution(t, prob, alg, nothing, nothing)
