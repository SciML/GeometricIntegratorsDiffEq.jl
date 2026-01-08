using GeometricIntegratorsDiffEq
using Test

using GeometricIntegrators
using DiffEqBase: solve, ReturnCode, SecondOrderODEProblem
import ODEProblemLibrary: prob_ode_2Dlinear

@testset "GeometricIntegratorsDiffEq" begin
    @testset "Standard ODE Problems" begin
        prob = prob_ode_2Dlinear

        sol = solve(prob, GIEuler(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIMidpoint(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIHeun2(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIHeun3(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRalston2(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRalston3(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRunge(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIKutta(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRK416(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRK438(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GISSPRK3(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GICrankNicolson(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIKraaijevangerSpijker(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIQinZhang(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GICrouzeix(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIImplicitEuler(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIImplicitMidpoint(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GISRK3(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIGLRK(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRadauIA(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GIRadauIIA(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIA(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIB(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIC(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIC̄(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIID(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIE(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIF(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "Second Order ODE Problems" begin
        # Second order ODE problem - use the new 5-argument convention with p parameter
        u0 = zeros(2)
        v0 = ones(2)
        f2 = function (dv, v, u, p, t)
            dv .= -u
        end
        prob = SecondOrderODEProblem{true}(f2, v0, u0, (0.0, 5.0))

        sol = solve(prob, GISymplecticEulerA(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GISymplecticEulerB(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIAIIIB(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        sol = solve(prob, GILobattoIIIBIIIA(2), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end
end # testset GeometricIntegratorsDiffEq

# Run ExplicitImports tests
if get(ENV, "GROUP", "All") == "All" || get(ENV, "GROUP", "") == "ExplicitImports"
    include("explicit_imports_test.jl")
end

# Run allocation tests if AllocCheck is available
if get(ENV, "GROUP", "All") == "All" || get(ENV, "GROUP", "") == "Nopre"
    include("alloc_tests.jl")
end

# Run JET static analysis tests
if get(ENV, "GROUP", "All") == "All" || get(ENV, "GROUP", "") == "JET"
    include("jet_tests.jl")
end
