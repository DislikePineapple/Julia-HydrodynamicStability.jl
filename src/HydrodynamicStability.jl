"""
    HydrodynamicStability

Developed by Sheng for solving eigenvalue problems for hydrodynamic stability systems, including LST, PSE and Receptivity problems.
It also need a systems for parameter research.
"""
module HydrodynamicStability

using ForwardDiff, LinearAlgebra, FFTW
using SpecialFunctions
using ProgressMeter

import Printf: @printf
import UnPack: @unpack
import RecipesBase: @recipe

# abstract type AbstractLinearProblem end
# abstract type AbstractNonlinearProblem end

## Problem
abstract type AbstractProblem end

abstract type AbstractPDEProblem <: AbstractProblem end
abstract type AbstractODEProblem <: AbstractProblem end
abstract type AbstractEVProblem <: AbstractProblem end
abstract type AbstractSCProblem <: AbstractProblem end

abstract type InstabilityProblem <: AbstractEVProblem end

## Algorithm
abstract type AbstractAlgorithm end
abstract type AbstractODEAlgorithm <: AbstractAlgorithm end
abstract type AbstractPDEAlgorithm <: AbstractAlgorithm end
abstract type AbstractEVPAlgorithm <: AbstractAlgorithm end

## Solution
abstract type AbstractSolution end

abstract type AbstractNonlinearSolution end
abstract type AbstractODESolution end
abstract type AbstractPDESolution end
abstract type AbstractEVPSolution end

function solve end
function initial end

struct NullParameter end

include("utils.jl")

include("ode/ode_problem.jl")
include("ode/ode_solution.jl")
include("ode/ode_algorithm.jl")
include("ode/ode_solve.jl")

include("nonlinear/nonlinear_problem.jl")
include("nonlinear/nonlinear_solution.jl")
include("nonlinear/nonlinear_solve.jl")

include("pde/pde_problem.jl")
include("pde/pde_solution.jl")
include("pde/pde_solve.jl")

export central_difference, simpsons_integral, chebyshev, chebyshevshift
export fft_expand, ifft_expand
export nest_vector, flatten_vector
export parameter_variation

export ODEProblem, BVProblem, NonlinearProblem
export ODESolution
export HeatProblem, NSHeatProblem
export solve

export Newton, Secant, Muller, Bisection, Falsi
export RK4, Shooting, FDM, SFDM, NSFDM
export PFDM, PSFDM, NPSFDM

end
