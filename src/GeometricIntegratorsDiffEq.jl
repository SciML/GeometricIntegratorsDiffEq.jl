module GeometricIntegratorsDiffEq

using Reexport
@reexport using DiffEqBase

using GeometricIntegrators

const warnkeywords =
    (:save_idxs, :d_discontinuities, :unstable_check, :save_everystep,
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

export GeometricIntegratorAlgorithm, GIEuler, GIMidpoint, GIHeun, GIKutta,
       GIKutta, GIERK4, GIERK438, GIImplicitEuler, GIImplicitMidpoint,
       GIRadIIA2, GIRadIIA3, GISRK3, GIGLRK, GISymplecticEulerA, GISymplecticEulerB,
       GILobattoIIIAIIIB2, GILobattoIIIBIIIA2

end # module
