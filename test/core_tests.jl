using GeometricIntegratorsDiffEq
using GeometricIntegrators
using DiffEqBase: solve, SecondOrderODEProblem
import ODEProblemLibrary: prob_ode_2Dlinear
using SciMLBase: successful_retcode
using Test

function test_successful_solve(prob, alg::GeometricIntegratorAlgorithm; dt)
    sol = solve(prob, alg; dt)
    @test successful_retcode(sol)
    return sol
end

@testset "Standard ODE Problems" begin
    prob = prob_ode_2Dlinear

    algorithms = GeometricIntegratorAlgorithm[
        GIEuler(), GIMidpoint(), GIHeun2(), GIHeun3(),
        GIRalston2(), GIRalston3(), GIRunge(), GIKutta(),
        GIRK4(), GIRK416(), GIRK438(), GISSPRK3(), GICrankNicolson(),
        GIKraaijevangerSpijker(), GIQinZhang(), GICrouzeix(),
        GIImplicitEuler(), GIImplicitMidpoint(), GISRK3(), GIGLRK(2),
        GIRadauIA(2), GIRadauIIA(2), GILobattoIIIA(2),
        GILobattoIIIB(2), GILobattoIIIC(2), GILobattoIIIC̄(2),
        GILobattoIIID(2), GILobattoIIIE(2), GILobattoIIIF(2),
    ]

    for alg in algorithms
        test_successful_solve(prob, alg; dt = 0.1)
    end
end

@testset "Second Order ODE Problems" begin
    # Second order ODE problem - use the new 5-argument convention with p parameter
    u0 = zeros(2)
    v0 = ones(2)
    f2 = function (dv, v, u, p, t)
        dv .= -u
    end
    prob = SecondOrderODEProblem{true}(f2, v0, u0, (0.0, 5.0))

    algorithms = GeometricIntegratorAlgorithm[
        GISymplecticEulerA(), GISymplecticEulerB(),
        GILobattoIIIAIIIB(2), GILobattoIIIBIIIA(2),
    ]

    for alg in algorithms
        test_successful_solve(prob, alg; dt = 0.1)
    end
end
