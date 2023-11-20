using HydrodynamicStability
using Aqua
# using Test

Aqua.test_all(HydrodynamicStability; ambiguities = false)

# @testset "Aqua tests (performance)" begin
#     # This tests that we don't accidentally run into
#     # https://github.com/JuliaLang/julia/issues/29393
#     Aqua.test_unbound_args(HydrodynamicStability)

#     # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
#     # Test that we're not introducing method ambiguities across deps
#     ambs = Aqua.detect_ambiguities(HydrodynamicStability; recursive = true)
#     pkg_match(pkgname, pkdir::Nothing) = false
#     pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
#     filter!(x -> pkg_match("HydrodynamicStability", pkgdir(last(x).module)), ambs)

#     # Uncomment for debugging:
#     # for method_ambiguity in ambs
#     #     @show method_ambiguity
#     # end
#     @warn "Number of method ambiguities: $(length(ambs))"
#     @test length(ambs) ≤ 13
# end

# @testset "Aqua tests (additional)" begin
#     Aqua.test_undefined_exports(HydrodynamicStability)
#     Aqua.test_stale_deps(HydrodynamicStability)
#     Aqua.test_deps_compat(HydrodynamicStability)
#     Aqua.test_project_extras(HydrodynamicStability)
#     # Aqua.test_project_toml_formatting(HydrodynamicStability)
#     # Aqua.test_piracy(HydrodynamicStability)
# end
