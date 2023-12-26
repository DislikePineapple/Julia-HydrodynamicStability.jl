abstract type ODEsAlgorithm <: AbstractAlgorithm end
struct RK4 <: ODEsAlgorithm end

abstract type BVProblemAlgorithm <: AbstractAlgorithm end
struct Shooting{I,IT} <: BVProblemAlgorithm
    ivp::I
    iter::IT
end

Shooting() = Shooting(RK4(), Secant())

struct FDM <: BVProblemAlgorithm end
struct FEM <: BVProblemAlgorithm end
