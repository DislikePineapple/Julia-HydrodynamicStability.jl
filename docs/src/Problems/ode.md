# ODE problems
Solve ordianary diffential equations with finit differential method.

## Boundary value problem
The general govening equation can be written as
```math
y' = f(y,t), 
```
which implies the equations are reduced to first-order vector form.
The above equations take the matrix form
```math
A\boldsymbol{y}'+ D\boldsymbol{y} =  F. 
```

## Finite Diffenrence Method
### Chebyshev method
### The discrete schemes
4-order center differential is used 
```math
  \frac{\mathrm{d}y_i}{\mathrm{d}t} = \frac{y_{i-2}-8y_{i-1}+8y_{i+1}-y_{i+2}}{12\Delta t}, 
```

### Example
Airy function
```math
y'' - t y = 0, 
```
which can be written as the first order form
```math
y_1' = y_2,
y_2' = ty_1.
```
Thesefore, in matrix form
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