# Backward finite differential method to solve parabilic PDE

## The eqution to be solved

The general form of parabilic equations follows:
$$
\begin{gather}
    V_{xx}\frac{\partial^2 u_{i,j}}{\partial x^2} + A\frac{\partial u_{i,j}}{\partial x} +\varGamma\frac{\partial u_{i,j}}{\partial t} + D u_{i,j} = F. 
\end{gather}
$$

## Algorithm, grid and discrete scheme: backward finite differential method
Temporal discretization: 2-order backward differential
$$ 
\begin{gather}
    \left(\frac{\partial u}{\partial t}\right)_{i,j} = \frac{3u_{i,j}-4u_{i-1,j}+u_{i-2,j}}{2\Delta t_{i,j}}, 
\end{gather}
$$
and reduce to the 1-order form at the entrance as
$$
\begin{gather}
    \left(\frac{\partial u}{\partial t}\right)_{i,j} = \frac{u_{i,j}-u_{i-1,j}}{\Delta t_{i,j}}. 
\end{gather}
$$
Spatial discretization: 5-point 4-order center differential
$$
\begin{gather}
\left(\frac{\partial u}{\partial x}\right)_{i,j} = \frac{u_{i,j-2}-8u_{i,j-1}+8u_{i,j+1}-u_{i,j+2}}{12\Delta x_j},\\
\left(\frac{\partial^2 u}{\partial x^2}\right)_{i,j} = \frac{-u_{i,j-2}+16u_{i,j-1}-30u_{i,j}+16u_{i,j+1}-u_{i,j+2}}{12\Delta x_j^2}. 
\end{gather} 
$$
The spatial discretization also need to reduce when close to the boundary.
<!-- ## Temp
The goved equation can be written as first order ODE:
$$ \frac{\partial u_{i,j}}{\partial x} = v_{i,j}, $$
$$ \frac{\partial v_{i,j}}{\partial x} = \frac{3\lambda }{2\Delta t} u_{i,j} -\frac{\lambda}{2\Delta t}(4u_{i-1,j}-u_{i-2,j}), $$
which can be written as the matrix form 
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
$$ -->

## Crank-Nicolson method

The mixed differential is used for the spatial discretization:
$$ 
\begin{gather}
\left(\frac{\partial u}{\partial t}\right)_{i,j} = \frac{u_{i,j}-u_{i-1,j}}{\Delta t_{i,j}},\\
\left(\frac{\partial u}{\partial x}\right)_{i,j} = \frac{u_{i,j+1}-u_{i,j-1}}{2\Delta x_{i,j}}+\frac{u_{i-1,j+1}-u_{i-1,j-1}}{2\Delta x_{i-1,j}}.
% \left(\frac{\partial^2 u}{\partial x^2}\right)_{i,j} = \frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{2\Delta x_{i,j}^2}+\frac{u_{i-1,j+1}-2u_{i-1,j}+u_{i-1,j-1}}{2\Delta x_{i-1,j}^2}.  
\end{gather}
$$
<!-- therefore the govening equation can be written as
$$ \frac{\partial^2 u_{i,j}}{\partial x^2} - \frac{3\lambda }{\Delta t} u_{i,j} = \frac{\lambda}{\Delta t}(4u_{i-1,j}-u_{i-2,j}) - \frac{\partial^2 u_{i-1,j}}{\partial x^2}, $$ -->

**Careful when deal with the boundary conditions.**

## Example
### Heat equation
The heat equation for the coefficient of heat conduction $\lambda=1$:
$$
\begin{gather}  
\frac{\partial^2 u}{\partial^2 x} - \lambda\frac{\partial u}{\partial t} = 0, 
\end{gather}
$$
Solve for $\lambda=1$ with the initial condition and boundary condition
$$
\begin{align}
    u(x,0) & = \sin^2(2\pi x)\quad \forall\quad 0\le x\le 1,\quad t\ge 0,\\ 
    u(0,t) & = 0 \qquad\qquad\ \ \forall\quad t>0,\\
    u(1,t) & = 0 \qquad\qquad\ \ \forall\quad t>0.
\end{align}
$$
The governing equations can be written in as the 1-order form
$$
\begin{gather}
    \frac{\partial u_1}{\partial x} - u_2 = 0,\\
    \frac{\partial u_2}{\partial x}  -\lambda\frac{\partial u_1}{\partial t}  = 0,
\end{gather}
$$
and in matrix form
$$
\begin{gather}
    A\frac{\partial u_{i,j}}{\partial x} +\varGamma\frac{\partial u_{i,j}}{\partial t} + D u_{i,j} = 0,
\end{gather}
$$
where
$$
\begin{gather}
    A=\begin{pmatrix}
        1 & 0 \\
        0 & 1
    \end{pmatrix},
    \quad
    \varGamma=\begin{pmatrix}
        0 & 0 \\
        -\lambda & 0
    \end{pmatrix},
    \quad
    D=\begin{pmatrix}
        0 & -1 \\
        0 & 0
    \end{pmatrix}.
\end{gather}
$$
