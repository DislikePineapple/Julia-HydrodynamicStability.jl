@doc """
Definde the ordinary differential equations solution

"""
struct ODESolution{uType,yType,P,A,uType2,RT}
    u::uType
    y::yType
    prob::P
    alg::A

    u_analytic::uType2
    retcode::RT
end

ODESolution(u, y, prob, alg) = ODESolution(u, y, prob, alg, nothing, nothing)
