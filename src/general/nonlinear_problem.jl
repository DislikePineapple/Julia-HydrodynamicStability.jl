struct NonlinaerProblem <: AbstractProblem
    f::F
    u0::uType
    p::P
    kwarg::K
end
