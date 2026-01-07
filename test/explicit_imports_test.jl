using ExplicitImports
using GeometricIntegratorsDiffEq
using Test

@testset "ExplicitImports" begin
    @testset "No implicit imports" begin
        @test check_no_implicit_imports(GeometricIntegratorsDiffEq) === nothing
    end

    @testset "No stale explicit imports" begin
        @test check_no_stale_explicit_imports(GeometricIntegratorsDiffEq) === nothing
    end
end
