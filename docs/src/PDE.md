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

### Streaky base flow
The linear equation for heating induced streaky baseflow:
$$ 
\begin{gather}
    2x \tilde{u}\frac{\partial \tilde{u}}{\partial x}-\left(2x\frac{\partial \tilde{f}}{\partial x}+\tilde{f}\right)\frac{\partial \tilde{u}}{\partial \eta}=\frac{\partial}{\partial \eta}\left(\tilde{C}\frac{\partial \tilde{u}}{\partial \eta}\right),\\\mathrm{with}\quad\frac{\partial \tilde{f}}{\partial \eta}=\tilde{u},
    \\
    \begin{split}    
        &2x \tilde{u}\frac{\partial \tilde{T}}{\partial x}-\left(2x\frac{\partial \tilde{f}}{\partial x}+\tilde{f}\right)\frac{\partial \tilde{T}}{\partial \eta}\\
    =   &\frac{1}{P_r}\frac{\partial}{\partial \eta}\left(\tilde{C}\frac{\partial \tilde{T}}{\partial \eta}\right)  +(\gamma-1)Ma^2\left(\frac{\partial \tilde{u}}{\partial \eta}\right)^2,
    \end{split}
    % \\
    % 2x\frac{\partial F_B}{\partial \eta}\frac{\partial T}{\partial x} - F_B \frac{\partial T}{\partial \eta} = \frac{1}{Pr}\frac{\partial^2 T}{\partial \eta^2},
\end{gather}
$$
with the initial and boundary conditions
$$
\begin{gather}
     \tilde{u}=1,\quad \tilde{T}=1\quad\mathrm{as}\quad\eta\to\infty,\\
    \tilde{f}=\tilde{u}=0,\quad \tilde{T}=T_w(x,\hat{z})\quad\mathrm{on}\quad\eta=0.
\end{gather}
$$
#### Linear equations when Chapman's law adopted
When Chapman's law is used, i.e. $C=1$, the governing equation can be written as 
$$
\begin{gather}
    f_B'''+f_Bf_B''=0,\\
    2xf_B'\frac{\partial \tilde{T}}{\partial x}-f_B\frac{\partial \tilde{T}}{\partial \eta}
    =   \frac{1}{P_r}\frac{\partial^2 \tilde{T}}{\partial \eta^2}  +(\gamma-1)Ma^2f_B''^2.
\end{gather}
$$
The D-H inverse transform can be solved together
$$
\begin{gather}
    \tilde{y}=\sqrt{2x}\int_{0}^{\eta}\tilde{T}(x,\eta^{\dag},\hat{z})\mathrm{d}\eta^{\dag} \quad\mathrm{with}\quad \tilde{y}=0\quad \mathrm{on}\quad\eta=0.
\end{gather}
$$
Furthemore, the energe equation follows the 1-order form
$$
\begin{align}
    \frac{\partial f_1}{\partial \eta} - f_2 & = 0 ,\\
    \frac{\partial f_2}{\partial \eta} - f_3 & = 0 ,\\
    \frac{\partial f_3}{\partial \eta}  & = F_3,\\
    \frac{\partial g_0}{\partial \eta} - \sqrt{2x} g_1 & = 0,\\
    \frac{\partial g_1}{\partial \eta} - g_2 & = 0,\\
    \frac{1}{P_r}\frac{\partial g_2}{\partial \eta} - 2x f_2\frac{\partial g_1}{\partial x}
    % -2x\frac{\partial f_1}{\partial x}g_2
    +f_1g_2 
    & = F_6,
\end{align}
$$
with
$$
\begin{align}
    F_3 & = -f_1f_3,\\
    F_6 & = -(\gamma-1)Ma^2f_3^2. 
\end{align}
$$
The initial and boundary conditions follows
$$
\begin{gather}
f_1\to0,\quad f_2\to0,\quad g_1\to T_B\quad\mathrm{as}\quad x\to0,\\
f_2\to1,\quad g_1\to1\quad\mathrm{as}\quad \eta\to\infty,\\
g_0=0,\quad g_1=T_w(x,z)\quad\mathrm{on}\quad \eta=0.
\end{gather}
$$
Appling the Fourier transform
$$
\begin{gather}
    \{f_1,f_2,f_3,F_3\}=\sum_{n=-\infty}^{\infty} \{\tilde{f}_{1,n},\tilde{f}_{2,n},\tilde{f}_{3,n},\tilde{F}_{3,n}\}\mathrm{e}^{\mathrm{i}n\beta \hat{z}},\\
    \{g_1,g_2,F_6\}=\sum_{n=-\infty}^{\infty} \{\tilde{g}_{1,n},\tilde{g}_{2,n},\tilde{F}_{2,n}\}\mathrm{e}^{\mathrm{i}n\beta \hat{z}}, \quad n\in \mathbb{Z},\\
    \mathrm{i.e.} \quad \tilde{F}_{3,n} = \mathscr{F}\left[-f_1f_3\right]_n,\quad\tilde{F}_{6,n} = \mathscr{F}\left[-(\gamma-1)Ma^2f_3^2\right]_n,
\end{gather}
$$
where $\mathscr{F}$ donate the Fourier operator.
Hence, we have
$$
\begin{align}
     \frac{\partial \tilde{f}_{1,n}}{\partial \eta} - \tilde{f}_{2,n} & = 0 ,\\
    \frac{\partial \tilde{f}_{2,n}}{\partial \eta} - \tilde{f}_{3,n} & = 0 ,\\
    \frac{\partial \tilde{f}_{3,n}}{\partial \eta}  & = \tilde{F}_{3,n},\\
    \frac{\partial \tilde{g}_{0,n}}{\partial \eta} - \sqrt{2x} \tilde{g}_{1,n} & = 0,\\
    \frac{\partial \tilde{g}_{1,n}}{\partial \eta} - \tilde{g}_{2,n} & =0 ,\\
    \frac{1}{P_r}\frac{\partial \tilde{g}_{2,n}}{\partial \eta} - 2x f_2\frac{\partial \tilde{g}_{1,n}}{\partial x}
    % -2x\frac{\partial f_1}{\partial x}g_2
    +f_1\tilde{g}_{2,n}
    & =
    \tilde{F}_{6,n}
    . 
\end{align}
$$
For the matrix form
$$
\begin{gather}
    A_n\frac{\partial u_{i,j}}{\partial \eta} +\varGamma_n\frac{\partial u_{i,j}}{\partial x} + D_n u_{i,j} = F_n,
\end{gather}
$$
where
$$
\begin{gather}
    A_n=\begin{bmatrix}
        1 & 0 & 0 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 & 0 & 0 \\
        0 & 0 & 0 & 0 & 1 & 0 \\
        0 & 0 & 0 & 0 & 0 & 1/Pr 
    \end{bmatrix},
    \\
    \varGamma_n=\begin{bmatrix}
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & -2xf_2 & 0 
    \end{bmatrix},
    \\
    D_n=\begin{bmatrix}
        0 & -1 & 0 & 0 & 0 & 0 \\
        0 & 0 & -1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 & -\sqrt{2x} & 0 \\
        0 & 0 & 0 & 0 & 0 & -1 \\
        0 & 0 & 0 & 0 & 0 & f_1 
    \end{bmatrix},
    \\
    \tilde{F}_n=\begin{bmatrix}
        0\\
        0\\
        -\mathscr{F}\left[f_1f_3\right]_n\\
        0\\
        0\\
        -\mathscr{F}\left[(\gamma-1)Ma^2f_3^2\right]_n
    \end{bmatrix}.
\end{gather}
$$
When the lower deck temperature is obtained, we need to carcualte the quantity in phsics spance by using the Dorodnitsyn-Howarth inverse transform
$$
\begin{gather}
    \tilde{y}=\sqrt{2x}\int_{0}^{\eta}\tilde{T}(x,\eta^{\dag},\hat{z})\mathrm{d}\eta^{\dag}.
\end{gather}
$$
The blowing velucity is carculated by
$$
\begin{gather}
    \tilde{v}=-\frac{1}{\sqrt{2x}}\left(f_B\tilde{T}-\int_0^\eta \tilde{T}\mathrm{d}\eta^{\dag}f_B'\right).
\end{gather}
$$

After that we need to solver the 