# GeometricIntegratorsDiffEq.jl

## Public API

```@docs
GeometricIntegratorAlgorithm
GIEuler
GIMidpoint
GIHeun2
GIHeun3
GIRalston2
GIRalston3
GIRunge
GIKutta
GIRK4
GIRK416
GIRK438
GISSPRK3
GICrankNicolson
GIKraaijevangerSpijker
GIQinZhang
GICrouzeix
GIImplicitEuler
GIImplicitMidpoint
GISRK3
GIGLRK
GIRadauIA
GIRadauIIA
GILobattoIIIA
GILobattoIIIB
GILobattoIIIC
GILobattoIIIC̄
GILobattoIIID
GILobattoIIIE
GILobattoIIIF
GISymplecticEulerA
GISymplecticEulerB
GILobattoIIIAIIIB
GILobattoIIIBIIIA
```

## Reexported SciML common interface

`using GeometricIntegratorsDiffEq` also brings in the parts of the SciML common interface
needed to build a problem, solve it with one of the methods above, and inspect the
result, so they do not have to be imported separately. These names are owned and
documented by [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/):

  - Problems: `ODEProblem`, `SecondOrderODEProblem`, `DynamicalODEProblem`,
    `EnsembleProblem`
  - Functions: `ODEFunction`, `DynamicalODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - `NullParameters`

The problem types are the ones these fixed-step methods accept: standard ODEs for the
Runge-Kutta wrappers, and second-order/dynamical problems for the symplectic and
partitioned wrappers. Callbacks and the iterator interface are deliberately *not*
reexported: `solve` errors when passed a `callback`, and this package implements no
`init`/`step!`. Anything else from SciMLBase must be imported from SciMLBase directly.
