module HydrodynamicStability

"""
## HydrodynamicStability

Developed by Sheng for solving eigenvalue problems for hydrodynamic stability systems, including LST, PSE and Receptivity problems.
It also need a systems for parameter research.
"""

# abstract type AbstractLinearProblem end
# abstract type AbstractNonlinearProblem end

abstract type AbstractBVProblem end
abstract type AbstractEVProblem end
abstract type AbstractSCProblem end

abstract type AbstractEVPSolution end

abstract type InstabilityProblem <: AbstractEVProblem end

end
