# Backward finite differential method to solve parabilic PDE

## The eqution to be solved
The heat equation:
$$ \frac{\partial^2 u}{\partial^2 x} - \lambda\frac{\partial u}{\partial t} = 0 . $$
The baseflow equations:
$$ 2x\frac{\partial F_B}{\partial \eta}\frac{\partial T}{\partial x} - F_B \frac{\partial T}{\partial \eta} = \frac{1}{Pr}\frac{\partial^2 T}{\partial \eta^2}. $$
The general equation:
$$ (Z_{xx}+\varGamma+A+D)u_{i,j}=F. $$

## Algorithm, grid and discrete scheme: backward finite differential method
Temporal discretization: 2-order backward differential
$$ \frac{\partial u}{\partial t} = \frac{3u_{i,j}-4u_{i-1,j}+u_{i-2,j}}{2\Delta t}. $$
Spatial discretization: 5-point 4-order center differential
$$ \frac{\partial u_{i,j}}{\partial x} = \frac{u_{i,j-2}-8u_{i,j-1}+8u_{i,j+1}-u_{i,j+2}}{12\Delta x}, $$
$$ \frac{\partial^2 u_{i,j}}{\partial x^2} = \frac{-u_{i,j-2}+16u_{i,j-1}-30u_{i,j}+16u_{i,j+1}-u_{i,j+2}}{12\Delta x^2}. $$
The goved equation can be written as first order ODE:
$$ \frac{\partial u_{i,j}}{\partial x} = v_{i,j}, $$
$$ \frac{\partial v_{i,j}}{\partial x} = \frac{3\lambda }{2\Delta t} u_{i,j} -\frac{\lambda}{2\Delta t}(4u_{i-1,j}-u_{i-2,j}), $$
which can be written as the matrix form
$$ A\frac{\partial u}{\partial x} + \varGamma\frac{\partial u}{\partial t} + D u+F=0,     $$
$$ I\boldsymbol{u}_{i,j}'=M\boldsymbol{u}_{i,j}+B, $$
where upscript $'$ is used to indicate the first order derivative, $\boldsymbol{u}_{i,j}=\{u_{i,j},v_{i,j}\}$ and
$$ M = \begin{bmatrix} 
0 & 1 \\
3\lambda/2\Delta t & 0
\end{bmatrix}, \quad B = \begin{bmatrix}
    0 \\ \lambda(4u_{i-1,j}-u_{i-2,j})/2\Delta t 
\end{bmatrix}. $$

The general governing equation takes the form
$$ V_{xx}\frac{\partial^2 u_{i,j}}{\partial x^2} + A\frac{\partial u_{i,j}}{\partial x} +\left(D+\frac{3\varGamma}{2\Delta t}\right) u_{i,j} = \frac{\varGamma}{2\Delta t}(4u_{i-1,j}-u_{i-2,j})+F. $$
The above equation could be written as the matrix form and the matrix takes
$$ M= \begin{bmatrix}
    0 & 1 \\
    -\frac{D}{V_{xx}}-\frac{3\varGamma}{2\Delta t V_{xx}}\ & - \frac{A}{V_{xx}}
\end{bmatrix}, \quad
B =\begin{bmatrix}
    0 \\
    \frac{\varGamma(4u_{i-1,j}-u_{i-2,j})}{2\Delta t V_{xx}}
\end{bmatrix}.
$$

## Crank-Nicolson method

The mixed differential is used for the spatial discretization:
$$ \frac{1}{2}\frac{\partial^2 u_{i,j}}{\partial x^2}+\frac{1}{2}\frac{\partial^2 u_{i-1,j}}{\partial x^2}, $$
therefore the govening equation can be written as
$$ \frac{\partial^2 u_{i,j}}{\partial x^2} - \frac{3\lambda }{\Delta t} u_{i,j} = \frac{\lambda}{\Delta t}(4u_{i-1,j}-u_{i-2,j}) - \frac{\partial^2 u_{i-1,j}}{\partial x^2}, $$

**Careful when deal with the boundary conditions.**