# Hydrodynamic Stability Toolkit
[![Build Status](https://github.com/DislikePineapple/Julia-HydrodynamicStability.jl/workflows/CI/badge.svg)](https://github.com/DislikePineapple/Julia-HydrodynamicStability.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/github/DislikePineapple/Julia-HydrodynamicStability.jl/branch/master/graph/badge.svg)](https://app.codecov.io/github/DislikePineapple/Julia-HydrodynamicStability.jl)
[![Global Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://dislikepineapple.github.io/Julia-HydrodynamicStability.jl/dev/)

## Overview

This package provides a comprehensive solution for eigenvalue problems associated with hydrodynamic stability systems.
It aims to establish a library encompassing Linear Stability Theory (LST), Parabolised Stability Equations (PSE), and Direct Numerical Systems (DNS).

## Structure

This package is designed based on the structure and fundamental concepts of [SciMLBase.jl](https://github.com/SciML/SciMLBase.jl), particularly in terms of the abstract types:

- Problem
- Algorithms
- Solution
- Function

Furthermore, the coding style adheres to the conventions established by SciML.
