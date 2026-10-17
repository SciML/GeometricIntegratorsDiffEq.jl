using GeometricIntegratorsDiffEq
using JET
using Test

@testset "JET static analysis" begin
    @testset "get_method_from_alg type stability" begin
        # Test that get_method_from_alg is type-stable for key algorithms
        # Using @test_opt to check for type instabilities

        # Explicit methods
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIEuler())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIMidpoint())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIRK4())

        # Implicit methods
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIImplicitEuler())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIImplicitMidpoint())

        # Parameterized methods
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIGLRK(2))
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GIRadauIA(2))
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GILobattoIIIA(2))

        # Symplectic methods
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GISymplecticEulerA())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.get_method_from_alg(GISymplecticEulerB())
    end

    @testset "requires_newton_solver type stability" begin
        # Test that requires_newton_solver is type-stable
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.requires_newton_solver(GIEuler())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.requires_newton_solver(GIMidpoint())
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.requires_newton_solver(GIRadauIA(2))
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIA(3))
        @test_opt target_modules = (GeometricIntegratorsDiffEq,) GeometricIntegratorsDiffEq.requires_newton_solver(GIImplicitMidpoint())
    end

    @testset "requires_newton_solver correctness" begin
        # Verify correct return values for requires_newton_solver
        # Methods that require Newton solver
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIRadauIA(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIRadauIIA(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIA(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIC(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIID(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIE(2)) == true
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIF(2)) == true

        # Methods that do not require Newton solver
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIEuler()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIMidpoint()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIRK4()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIImplicitEuler()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GIImplicitMidpoint()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GISymplecticEulerA()) == false
        @test GeometricIntegratorsDiffEq.requires_newton_solver(GISymplecticEulerB()) == false
    end
end
