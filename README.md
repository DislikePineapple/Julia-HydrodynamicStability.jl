# Hydrodynamic Stability Toolkit
[![Global Docs](https://img.shields.io/badge/docs-HST-blue.svg)](https://dislikepineapple.github.io/Julia-HydrodynamicStability.jl/stable/)

## Introduction

This is a package for solving eigenvalue problems of hydrodynamics stability systems. 
This package aims to establish a lib that includes linear stability theory (LST), partial stability equations (PSE) and direct numerical system (DNS).

## Abstract type

The present package uses the structure of [SciMLBase.jl](https://github.com/SciML/SciMLBase.jl) for the abstract type:

- Problem
- Algorithms
- Solution
- Function

Also, the code style follows SciML.