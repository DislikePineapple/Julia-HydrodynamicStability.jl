@doc """
Definde an nonlinear eigenvalue problem, which is used for the instabiliry analysis.
L(α) = 0,
where L is the operator, e.g. O-S operator, and α is the eigenvalue for a spatial problem.

Especially for my problem, it can't be write as a generalized eigenvalue problem.

    f: the operator to be solved
    v0: the guess of eigenvalue for the problem
"""
struct EVProblem{F,V,P,K} <: AbstractEVProblem
    f::F
    v0::V

    p::P
    kwarg::K

    function EVProblem(f, v0, p = NullParameters; kwarg...)
        new{typeof(f),typeof(v0),typeof(p),typeof(kwarg)}(f, v0, p, kwarg)
    end
end
