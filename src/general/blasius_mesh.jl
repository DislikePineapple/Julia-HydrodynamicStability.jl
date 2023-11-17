struct BlasiusMesh{Coordinate<:AbstractArray}
    η::Coordinate
    # y::Coordinate

    Nη::Int
    η∞::Any

    function BlasiusMesh(Nη, η∞)
        η = range(0, η∞, Nη)
        new{typeof(η)}(η, Nη, η∞)
    end
end
