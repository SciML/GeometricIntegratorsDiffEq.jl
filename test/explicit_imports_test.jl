using ExplicitImports
using GeometricIntegratorsDiffEq
using Test

@testset "Explicit Imports" begin
    @test check_no_implicit_imports(GeometricIntegratorsDiffEq) === nothing
    @test check_no_stale_explicit_imports(GeometricIntegratorsDiffEq) === nothing
end
