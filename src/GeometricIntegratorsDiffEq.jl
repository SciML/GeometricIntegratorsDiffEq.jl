module GeometricIntegratorsDiffEq

using Reexport: Reexport, @reexport
@reexport using DiffEqBase: DiffEqBase
using SciMLBase: SciMLBase, ReturnCode, check_keywords, isinplace, warn_compat

using GeometricIntegrators: GeometricIntegrators, CrankNicolson, Crouzeix,
    ExplicitEuler, ExplicitMidpoint, Gauss, Heun2, Heun3, ImplicitEuler,
    ImplicitMidpoint, KraaijevangerSpijker, Kutta3, LobattoIIIA,
    LobattoIIIAIIIB, LobattoIIIB, LobattoIIIBIIIA, LobattoIIIC, LobattoIIID,
    LobattoIIIE, LobattoIIIF, Newton, QinZhang, RK4, RK416, RK438, RadauIA,
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

include("algorithms.jl")
include("solve.jl")

export GeometricIntegratorAlgorithm, GIEuler, GIMidpoint, GIHeun2, GIHeun3,
    GIRalston2, GIRalston3, GIRunge, GIKutta, GIRK4, GIRK416, GIRK438, GISSPRK3,
    GICrankNicolson, GIKraaijevangerSpijker, GIQinZhang, GICrouzeix,
    GIImplicitEuler, GIImplicitMidpoint, GISRK3,
    GIGLRK, GILobattoIIIA, GILobattoIIIB, GILobattoIIIC, GILobattoIIIC̄,
    GILobattoIIID, GILobattoIIIE, GILobattoIIIF, GIRadauIA, GIRadauIIA,
    GISymplecticEulerA, GISymplecticEulerB, GILobattoIIIAIIIB,
    GILobattoIIIBIIIA

end # module
