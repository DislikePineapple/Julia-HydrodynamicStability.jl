struct FreeStream
    Ma::Float64
    Te::Float64
    Re::Float64
end

struct Grid
    x::AbstractArray
    y::AbstractArray
    z::AbstractArray
end

mutable struct Wave
    ω::Complex
    α::Complex
    β::Complex
end

mutable struct Heating{Ffunction<:Function,Sfunction<:Function}
    h::Float64
    β::Float64
    f::Ffunction
    s::Sfunction
    function Heating(h, β)
        d = 0.3
        x_0 = 0.955
        f(x) = x * exp(-(x - x_0)^2 / d^2)
        s(z) = cos(β * z)
        new{typeof(f),typeof(s)}(h, β, f, s)
    end
end
