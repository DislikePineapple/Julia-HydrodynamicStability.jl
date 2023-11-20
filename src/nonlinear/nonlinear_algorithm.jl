abstract type NonlinearAlgorithm <: AbstractAlgorithm end
struct Bisection <: NonlinearAlgorithm end
struct Falsi <: NonlinearAlgorithm end
struct Muller <: NonlinearAlgorithm end
struct Secant <: NonlinearAlgorithm end
