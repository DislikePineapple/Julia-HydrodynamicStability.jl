# Nonlinear Equations Problem

## Nonlinear Problem

```@docs
NonlinearProblem{F, DF, tType, P, K}
```

## Algorithm

### Bisection method

if ``f(x)`` is a continuity function in range ``[a, b]``, and which satisfy ``f(a)f(b)<0``, then there exist a roof between ``a`` and ``b``, which means there is a ``x_0`` satisfing ``a<x_0<b`` and ``f(x_0)=0``.

Then we can check the value of ``f(c)``, where ``c=\frac{a+b}{2}``, and determine whether ``f(a)f(c)<0`` or ``f(c)f(b)<0``, untill the value of ``f(c)<\mathrm{tol}``.

### Newton-Raphson method

Newton-Raphson method is a root finding method for the nonlinear equations (1). First, we give a initial guess ``x_0`` and we have
```math
    \begin{gather}
        f'(x_0)=\frac{f(x_0)-0}{x_0-x_1},\\
        x_1=x_0-\frac{f(x_0)}{f'(x_0)},
    \end{gather}
```
where ``x_1`` is the approximate root we are looking for. Repeat the progress of (3) untill the solution converges, we have the iteration equation that
```math
    x_0 = \mathrm{intial\ guess},\\
    x_{i+1}=x_i-\frac{f(x_i)}{f'(x_i)},\quad i=0,1,2,\cdots
```

### Muller's method