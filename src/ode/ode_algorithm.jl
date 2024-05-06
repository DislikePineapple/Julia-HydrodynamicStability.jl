abstract type ODEAlgorithm <: AbstractODEAlgorithm end
struct RK4 <: ODEAlgorithm end

abstract type BVPAlgorithm <: AbstractODEAlgorithm end
struct Shooting{I, IT} <: BVPAlgorithm
    ivp::I
    iter::IT
end

Shooting() = Shooting(RK4(), Secant())

struct FDM <: BVPAlgorithm end
struct SFDM <: BVPAlgorithm end
struct NSFDM <: BVPAlgorithm end

struct FEM <: BVPAlgorithm end
