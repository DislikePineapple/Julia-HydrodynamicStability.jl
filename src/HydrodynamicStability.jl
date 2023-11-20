module HydrodynamicStability

"""
## HydrodynamicStability

Developed by Sheng for solving eigenvalue problems for hydrodynamic stability systems, including LST, PSE and Receptivity problems.
It also need a systems for parameter research.
"""

import UnPack: @unpack
import RecipesBase: @recipe

# abstract type AbstractLinearProblem end
# abstract type AbstractNonlinearProblem end

abstract type AbstractProblem end

abstract type AbstractODEProblem <: AbstractProblem end
abstract type AbstractBVProblem <: AbstractProblem end
abstract type AbstractEVProblem <: AbstractProblem end
abstract type AbstractSCProblem <: AbstractProblem end

abstract type InstabilityProblem <: AbstractEVProblem end

abstract type AbstractAlgorithm end

abstract type AbstractODESolution end

function solve end
function initial end

include("utils.jl")

include("general/ode_problem.jl")
include("general/ode_solution.jl")
include("general/ode_algorithm.jl")
include("general/solve.jl")

export ODEProblem
export solve

export RK4

end
