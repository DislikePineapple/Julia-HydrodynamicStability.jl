# Ordinary Differential Equations
Solve ordinary differential equations with finite differential method.

## Problems

### Initial value problem

```@docs
ODEProblem{F, uType, yType, P, K}
```

### Boundary value problem

The general governing equation can be written as

```math
    u' = f(u,t), 
```

with the boundary conditions

```math
    u(0) = a,\quad u(1) = b.
```

Following the mathematical definition, the numerical specification is given as

```@docs
BVProblem{F, BC, U0, Y, P, K}
```

#### Linear boundary value problem

For the linear problem,
the equations can be reduced to first-order vector form, which take the matrix form

```math
    A\boldsymbol{u(t)}' + D\boldsymbol{u(t)} = F(t). 
```

#### Nonlinear boundary value problem
For the nonlinear problem, the equations can be reduced to first-order vector form, which take the matrix form

```math
    A\boldsymbol{u}(t)' + D\boldsymbol{u}(t) = N(\boldsymbol{u},t) + F(t). 
```
<!-- Here, the superscript and subscript choose as a n*n Matrix with subscript i row, j column, and superscript k donate the iteration steps. -->
Then, in order to use the New-Raphson method to solve the nonlinear problem ``R(\boldsymbol{u}^k)=0``, the equations should be written as the Jacobian matrix that

```math
    \boldsymbol{u}^{k+1} = \boldsymbol{u}^k-J^{-1}(\boldsymbol{u}^k)R(\boldsymbol{u}^k),
```

where

```math
    R \begin{bmatrix}\boldsymbol{u}_1 \\ \boldsymbol{u}_2 \\ \vdots \\ \boldsymbol{u}_n \end{bmatrix} = M\boldsymbol{u}-\begin{bmatrix}N(\boldsymbol{u}_1,t_1) + F(t_1) \\ N(\boldsymbol{u}_2,t_2) + F(t_2) \\ \vdots \\ N(\boldsymbol{u}_n,t_n) + F(t_n)\end{bmatrix},

```

with

```math
    M = \mathrm{kron}(Diff,A)+\mathrm{kron}(I,D),
```

and

```math
    J(\boldsymbol{u}) = M-\begin{bmatrix}N'(\boldsymbol{u}_1,t_1) & 0 & \cdots & 0  \\ 0 & N'(\boldsymbol{u}_2,t_2) & \cdots & 0  \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & N'(\boldsymbol{u}_n,t_n) \end{bmatrix}.
```

Then, to deal with the boundary conditions,


The shooting method can be used to solve the nonlinear boundary value problem as well.

## Solution 

```@docs
ODESolution{uType, yType, P, A, uType2, RT}
```

## Algorithm

### 4th-order Runge-Kutta method

### Shooting Method

### Finite Differential Method

#### Chebyshev collection

#### The discrete schemes
4-order center differential is used 
```math
  \frac{\mathrm{d}y_i}{\mathrm{d}t} = \frac{y_{i-2}-8y_{i-1}+8y_{i+1}-y_{i+2}}{12\Delta t}, 
```

## Example

### Airy function

```math
u'' - t u = 0,
```
with the boundary condition
```math
u(0) = \operatorname{Ai}(0),\\
u(10) = \operatorname{Ai}(10),
```
where ``\operatorname{Ai}(x)`` is the Airy function of the first kind.

The equation can be written as the first order form
```math
    y_1' = y_2,\\
    y_2' = ty_1.
```
Furthmore, write it in the matrix form
```math
    \varGamma\frac{\mathrm{d} y_{i}}{\mathrm{d} t} + D y_{i} = 0,
```
where
```math
    \varGamma=\begin{pmatrix}
        1 & 0 \\
        0 & 1
    \end{pmatrix},
    \quad
    D=\begin{pmatrix}
        0  & -1 \\
        -t & 0
    \end{pmatrix}.
```

### Nonlinear boundary value problem

#### Mathematical Specification of the Nonlinear Problem

```math
    u'' = u - u^2,
```

with the boundary condition

```math
    u(0) = 1,\quad u(1) = 4.
```

#### Numerical Specification of the Nonlinear Problem

The equation can be write in the matrix form

```math
    A\boldsymbol{u}(t)' + D\boldsymbol{u}(t) = N(\boldsymbol{u},t). 
```

where the coefficient matrix are

```math
    A = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix},\quad D = \begin{bmatrix} 0 & -1 \\ -1 & 0 \end{bmatrix},
```

and the nonlinear term is

```math
    N(\boldsymbol{u},t) = \begin{bmatrix} 0 \\ - u^2 \end{bmatrix}, \quad 
    DN(\boldsymbol{u},t) = \begin{bmatrix} 0 & 0  \\ -2u & 0 \end{bmatrix}.
```

### Blasius equation

#### Mathematical Specification of the Blasius Equation
solve the Blasius equation

```math
    f''' + \frac{1}{2} f f'' = 0,
```

with the boundary condition

```math
    f(0) = f'(0) = 0, \quad 
    f' \to 1 \quad \mathrm{as} \quad y \to \infty.
```

#### Numerical Specification of the Blasius Equation

**Numerical algorithm: Shooting method**

The equation can be write as the first order ODE system

```julia
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -1 / 2 * u[1] * u[3]
```

with the boundary condition

```julia
    residual[1] = u[begin][1]
    residual[2] = u[begin][2]
    residual[3] = u[end][2] - 1
```

The domain is `[0, 15]`, and the initial guess is `[0, 0, 0.3]`.

**Numerical algorithm: NSFDM**

The equations can be write in the matrix form

```math
    A\boldsymbol{u}(t)' + D\boldsymbol{u}(t) = N(\boldsymbol{u},t). 
```

where the coefficient matrix are

```math
    A = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix},\quad 
    D = \begin{bmatrix} 0 & -1 & 0 \\ 0 & 0 & -1 \\ 0 & 0 & 0 \end{bmatrix},
```

and the nonlinear term is

```math
    N(\boldsymbol{u},t) = -\frac{1}{2} \begin{bmatrix} 0 \\ 0 \\ f f'' \end{bmatrix}, \quad 
    DN(\boldsymbol{u},t) = -\frac{1}{2} \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ f'' & 0 & f \end{bmatrix}.
```