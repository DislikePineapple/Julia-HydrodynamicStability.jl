struct Shooting{I,IV}
    ivp::I
    iter::IT
end

Shooting() = Shooting(RK4(), Muller())

function SciMLBase.solve(prob::BVProblem, alg, args...; kwargs...)
    solve(init(prob, alg, args...; kwargs...))
end

function init(prob::BVProblem, alg, args...; kwargs...)

end
