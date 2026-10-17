module GeometricIntegratorsDiffEq

using DiffEqBase: DiffEqBase
using SciMLBase: SciMLBase
using SciMLLogging: @SciMLMessage
using RecursiveArrayTools: RecursiveArrayTools
using SimpleSolvers: SimpleSolvers

# The SciML common interface that GeometricIntegratorsDiffEq reexports (see the second
# `export` below), so that `using GeometricIntegratorsDiffEq` on its own is enough to
# build a problem, solve it with one of the GeometricIntegrators methods, and inspect
# the result -- which is exactly what the README and the docstring examples do. The
# problem types are the ones these fixed-step methods accept: standard ODEs and the
# dynamical/second-order problems the symplectic and partitioned methods integrate.
# Callbacks and the iterator interface are deliberately absent: `solve` errors on a
# `callback`, and this package implements no `init`/`step!`. Every name stays owned and
# documented upstream.
using SciMLBase: DynamicalODEFunction, DynamicalODEProblem, EnsembleAnalysis,
    EnsembleDistributed, EnsembleProblem, EnsembleSerial, EnsembleSolution,
    EnsembleSplitThreads, EnsembleSummary, EnsembleThreads, NullParameters, ODEFunction,
    ODEProblem, ODESolution, ReturnCode, SecondOrderODEProblem, remake, solve,
    successful_retcode

using GeometricIntegrators: GeometricIntegrators, CrankNicolson, Crouzeix,
    ExplicitEuler, ExplicitMidpoint, Gauss, Heun2, Heun3, ImplicitEuler,
    ImplicitMidpoint, KraaijevangerSpijker, Kutta3, LobattoIIIA,
    LobattoIIIAIIIB, LobattoIIIB, LobattoIIIBIIIA, LobattoIIIC, LobattoIIID,
    LobattoIIIE, LobattoIIIF, QinZhang, RK4, RK416, RK438, RadauIA,
    RadauIIA, Ralston2, Ralston3, Runge2, SRK3, SSPRK3, SymplecticEulerA,
    SymplecticEulerB, integrate

const warnkeywords = (
    :save_idxs, :d_discontinuities, :unstable_check, :save_everystep,
    :save_end, :initialize_save, :adaptive, :abstol, :reltol, :dtmax,
    :dtmin, :force_dtmin, :internalnorm, :gamma, :beta1, :beta2,
    :qmax, :qmin, :qsteady_min, :qsteady_max, :qoldinit, :failfactor,
    :maxiters, :isoutofdomain, :unstable_check,
    :calck, :progress, :tstops, :saveat, :dense,
)

function __init__()
    return global warnlist = Set(warnkeywords)
end

newton_solver_method() = if isdefined(SimpleSolvers, :NewtonMethod)
    SimpleSolvers.NewtonMethod()
else
    SimpleSolvers.Newton()
end

include("algorithms.jl")
include("solve.jl")

using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    @compile_workload begin
        GIEuler()
        GIMidpoint()
        GIHeun2()
        GIGLRK(2)
        newton_solver_method()
    end
end

export GeometricIntegratorAlgorithm, GIEuler, GIMidpoint, GIHeun2, GIHeun3,
    GIRalston2, GIRalston3, GIRunge, GIKutta, GIRK4, GIRK416, GIRK438, GISSPRK3,
    GICrankNicolson, GIKraaijevangerSpijker, GIQinZhang, GICrouzeix,
    GIImplicitEuler, GIImplicitMidpoint, GISRK3,
    GIGLRK, GILobattoIIIA, GILobattoIIIB, GILobattoIIIC, GILobattoIIIC̄,
    GILobattoIIID, GILobattoIIIE, GILobattoIIIF, GIRadauIA, GIRadauIIA,
    GISymplecticEulerA, GISymplecticEulerB, GILobattoIIIAIIIB,
    GILobattoIIIBIIIA

# Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
export DynamicalODEFunction, DynamicalODEProblem, EnsembleAnalysis, EnsembleDistributed,
    EnsembleProblem, EnsembleSerial, EnsembleSolution, EnsembleSplitThreads,
    EnsembleSummary, EnsembleThreads, NullParameters, ODEFunction, ODEProblem,
    ODESolution, ReturnCode, SecondOrderODEProblem, remake, solve, successful_retcode

end # module
