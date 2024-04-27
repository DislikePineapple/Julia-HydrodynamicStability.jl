struct NonlinearProblem{F, tType, P, K} <: AbstractProblem
    f::F
    t0::tType

    p::P
    kwarg::K

    function NonlinearProblem(f, t0, p = NullParameter(); kwargs...)
        new{typeof(f), typeof(t0), typeof(p), typeof(kwargs)}(f, t0, p, kwargs)
    end
end
