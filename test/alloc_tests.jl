using GeometricIntegratorsDiffEq
using AllocCheck
using Test

@testset "Allocation-free helper functions" begin
    # Test that algorithm method dispatch is allocation-free
    @testset "get_method_from_alg" begin
        @test begin
            @check_allocs checkf1(::GIEuler) = GeometricIntegratorsDiffEq.get_method_from_alg(GIEuler())
            true
        end

        @test begin
            @check_allocs checkf2(::GIMidpoint) = GeometricIntegratorsDiffEq.get_method_from_alg(GIMidpoint())
            true
        end

        @test begin
            @check_allocs checkf3(::GIImplicitMidpoint) = GeometricIntegratorsDiffEq.get_method_from_alg(GIImplicitMidpoint())
            true
        end

        @test begin
            @check_allocs checkf4(::GIRK4) = GeometricIntegratorsDiffEq.get_method_from_alg(GIRK4())
            true
        end

        @test begin
            @check_allocs checkf5(::GISymplecticEulerA) = GeometricIntegratorsDiffEq.get_method_from_alg(GISymplecticEulerA())
            true
        end
    end

    # Test that Newton solver check is allocation-free
    @testset "requires_newton_solver" begin
        @test begin
            @check_allocs checkns1(::GIEuler) = GeometricIntegratorsDiffEq.requires_newton_solver(GIEuler())
            true
        end

        @test begin
            @check_allocs checkns2(::GIRadauIA) = GeometricIntegratorsDiffEq.requires_newton_solver(GIRadauIA(2))
            true
        end

        @test begin
            @check_allocs checkns3(::GILobattoIIIA) = GeometricIntegratorsDiffEq.requires_newton_solver(GILobattoIIIA(3))
            true
        end

        @test begin
            @check_allocs checkns4(::GIImplicitMidpoint) = GeometricIntegratorsDiffEq.requires_newton_solver(GIImplicitMidpoint())
            true
        end
    end
end
