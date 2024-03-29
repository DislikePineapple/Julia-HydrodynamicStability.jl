"""
    HydrodynamicStability

Developed by Sheng for solving eigenvalue problems for hydrodynamic stability systems, including LST, PSE and Receptivity problems.
It also need a systems for parameter research.
"""
module HydrodynamicStability

using ForwardDiff, LinearAlgebra
# using FiniteDiff

import Printf: @printf
import UnPack: @unpack
import RecipesBase: @recipe

# abstract type AbstractLinearProblem end
# abstract type AbstractNonlinearProblem end

## Problem
abstract type AbstractProblem end

abstract type AbstractODEProblem <: AbstractProblem end
abstract type AbstractBVProblem <: AbstractProblem end
abstract type AbstractEVProblem <: AbstractProblem end
abstract type AbstractSCProblem <: AbstractProblem end

abstract type InstabilityProblem <: AbstractEVProblem end

## Algorithm
abstract type AbstractAlgorithm end
abstract type AbstractEVPAlgorithm <: AbstractAlgorithm end

## Solution
abstract type AbstractSolution end

abstract type AbstractNonlinearSolution end
abstract type AbstractODESolution end
abstract type AbstractEVPSolution end

function solve end
function initial end

struct NullParameter end

include("ode/ode_problem.jl")
include("ode/ode_solution.jl")
include("ode/ode_algorithm.jl")
include("ode/ode_solve.jl")

include("nonlinear/nonlinear_problem.jl")
include("nonlinear/nonlinear_solution.jl")
include("nonlinear/nonlinear_solve.jl")

export ODEProblem, BVProblem, NonlinearProblem
export solve

export Secant, Muller, Bisection, Falsi
export RK4, Shooting, FDM

end
