module GeometricIntegratorsDiffEq

using Reexport: @reexport
@reexport using DiffEqBase

using DiffEqBase: check_keywords, warn_compat
using GeometricIntegrators

const warnkeywords = (:save_idxs, :d_discontinuities, :unstable_check, :save_everystep,
    :save_end, :initialize_save, :adaptive, :abstol, :reltol, :dtmax,
    :dtmin, :force_dtmin, :internalnorm, :gamma, :beta1, :beta2,
    :qmax, :qmin, :qsteady_min, :qsteady_max, :qoldinit, :failfactor,
    :maxiters, :isoutofdomain, :unstable_check,
    :calck, :progress, :tstops, :saveat, :dense)

function __init__()
    global warnlist = Set(warnkeywords)
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
