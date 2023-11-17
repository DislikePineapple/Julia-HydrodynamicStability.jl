@doc """
Definde the ordinary differential equations solution

"""
struct ODESolution{uType,tType,P,A,uType2}
    u::uType
    t::tType
    prob::P
    alg::A

    u_analytic::uType2
    retcode::ReturnCode.T
end

ODESolution(u, t, prob, alg) = ODESolution(u, t, prob, alg, nothing, nothing)
