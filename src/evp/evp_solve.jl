struct NonEVP <: AbstractEVPAlgorithm end

function solve(prob::EVProblem, alg::NonEVP, arg...; kwarg...) end
