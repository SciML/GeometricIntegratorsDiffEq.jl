"""
    GeometricIntegratorAlgorithm <: SciMLBase.AbstractODEAlgorithm

Abstract interface for the fixed-step GeometricIntegrators.jl methods exposed through
the SciML `solve` interface.

Concrete subtypes are algorithm choices passed as the second argument to `solve`. They
wrap tableaus and methods from GeometricIntegrators.jl while returning a SciML solution
object.

# Interface

- `solve(prob, alg::GeometricIntegratorAlgorithm; dt, kwargs...)` requires a fixed
    timestep `dt`; these algorithms do not choose timesteps adaptively.
- Standard ODE algorithms support `SciMLBase.AbstractODEProblem`s whose problem type is
    `SciMLBase.StandardODEProblem`.
- Symplectic and partitioned algorithms support second-order or dynamical ODE problems
    through `SciMLBase.AbstractDynamicalODEProblem`.
- Callback handling is not supported. Passing `callback` throws an error.
- Common adaptive or dense-output keywords that GeometricIntegrators.jl cannot honor are
    ignored with a SciML logging warning controlled by `verbose`.
- The public extension surface is the `SciMLBase.AbstractODEAlgorithm` contract. Adding
    new wrappers currently also requires package-internal method mappings, so downstream
    packages should not subtype `GeometricIntegratorAlgorithm` without coordinating that
    mapping in GeometricIntegratorsDiffEq.

# Keywords

- `dt::Number`: Required fixed timestep used by the wrapped GeometricIntegrators.jl
    method.
- `save_start::Bool = true`: Include the initial condition in the returned solution.
- `timeseries_errors::Bool = true`: Forwarded to `SciMLBase.build_solution`.
- `alias_u0::Bool = false`: Reuse `prob.u0` for mutable array initial conditions instead
    of copying it before calling the wrapped integrator.
- `verbose = true`: Controls warnings for unsupported solver keywords.

# Examples

```julia
using DiffEqBase, GeometricIntegratorsDiffEq

f(u, p, t) = -u
prob = ODEProblem(f, 1.0, (0.0, 1.0))
sol = solve(prob, GIEuler(); dt = 0.1)
```

```julia
using DiffEqBase, GeometricIntegratorsDiffEq

function acceleration!(dv, v, u, p, t)
    dv .= -u
end

prob = SecondOrderODEProblem{true}(acceleration!, ones(2), zeros(2), (0.0, 1.0))
sol = solve(prob, GISymplecticEulerA(); dt = 0.1)
```
"""
abstract type GeometricIntegratorAlgorithm <: SciMLBase.AbstractODEAlgorithm end

"""
    GIEuler()

Fixed-step explicit Euler wrapper for `GeometricIntegrators.ExplicitEuler`.

# Fields

This type has no public fields.

# Usage

Pass `GIEuler()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIEuler <: GeometricIntegratorAlgorithm end

"""
    GIMidpoint()

Fixed-step explicit midpoint wrapper for `GeometricIntegrators.ExplicitMidpoint`.

# Fields

This type has no public fields.

# Usage

Pass `GIMidpoint()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIMidpoint <: GeometricIntegratorAlgorithm end

"""
    GIHeun2()

Fixed-step Heun second-order wrapper for `GeometricIntegrators.Heun2`.

# Fields

This type has no public fields.

# Usage

Pass `GIHeun2()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIHeun2 <: GeometricIntegratorAlgorithm end

"""
    GIHeun3()

Fixed-step Heun third-order wrapper for `GeometricIntegrators.Heun3`.

# Fields

This type has no public fields.

# Usage

Pass `GIHeun3()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIHeun3 <: GeometricIntegratorAlgorithm end

"""
    GIRalston2()

Fixed-step Ralston second-order wrapper for `GeometricIntegrators.Ralston2`.

# Fields

This type has no public fields.

# Usage

Pass `GIRalston2()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRalston2 <: GeometricIntegratorAlgorithm end

"""
    GIRalston3()

Fixed-step Ralston third-order wrapper for `GeometricIntegrators.Ralston3`.

# Fields

This type has no public fields.

# Usage

Pass `GIRalston3()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRalston3 <: GeometricIntegratorAlgorithm end

"""
    GIRunge()

Fixed-step Runge second-order wrapper for `GeometricIntegrators.Runge2`.

# Fields

This type has no public fields.

# Usage

Pass `GIRunge()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRunge <: GeometricIntegratorAlgorithm end

"""
    GIKutta()

Fixed-step Kutta third-order wrapper for `GeometricIntegrators.Kutta3`.

# Fields

This type has no public fields.

# Usage

Pass `GIKutta()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIKutta <: GeometricIntegratorAlgorithm end

"""
    GIRK4()

Fixed-step classical fourth-order Runge-Kutta wrapper for `GeometricIntegrators.RK4`.

# Fields

This type has no public fields.

# Usage

Pass `GIRK4()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRK4 <: GeometricIntegratorAlgorithm end

"""
    GIRK416()

Fixed-step RK416 wrapper for `GeometricIntegrators.RK416`.

# Fields

This type has no public fields.

# Usage

Pass `GIRK416()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRK416 <: GeometricIntegratorAlgorithm end

"""
    GIRK438()

Fixed-step RK438 wrapper for `GeometricIntegrators.RK438`.

# Fields

This type has no public fields.

# Usage

Pass `GIRK438()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRK438 <: GeometricIntegratorAlgorithm end

"""
    GISSPRK3()

Fixed-step third-order strong-stability-preserving Runge-Kutta wrapper for
`GeometricIntegrators.SSPRK3`.

# Fields

This type has no public fields.

# Usage

Pass `GISSPRK3()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GISSPRK3 <: GeometricIntegratorAlgorithm end

"""
    GICrankNicolson()

Fixed-step Crank-Nicolson wrapper for `GeometricIntegrators.CrankNicolson`.

# Fields

This type has no public fields.

# Usage

Pass `GICrankNicolson()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GICrankNicolson <: GeometricIntegratorAlgorithm end

"""
    GIKraaijevangerSpijker()

Fixed-step Kraaijevanger-Spijker wrapper for
`GeometricIntegrators.KraaijevangerSpijker`.

# Fields

This type has no public fields.

# Usage

Pass `GIKraaijevangerSpijker()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIKraaijevangerSpijker <: GeometricIntegratorAlgorithm end

"""
    GIQinZhang()

Fixed-step Qin-Zhang wrapper for `GeometricIntegrators.QinZhang`.

# Fields

This type has no public fields.

# Usage

Pass `GIQinZhang()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIQinZhang <: GeometricIntegratorAlgorithm end

"""
    GICrouzeix()

Fixed-step Crouzeix wrapper for `GeometricIntegrators.Crouzeix`.

# Fields

This type has no public fields.

# Usage

Pass `GICrouzeix()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GICrouzeix <: GeometricIntegratorAlgorithm end

"""
    GIImplicitEuler()

Fixed-step implicit Euler wrapper for `GeometricIntegrators.ImplicitEuler`.

# Fields

This type has no public fields.

# Usage

Pass `GIImplicitEuler()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIImplicitEuler <: GeometricIntegratorAlgorithm end

"""
    GIImplicitMidpoint()

Fixed-step implicit midpoint wrapper for `GeometricIntegrators.ImplicitMidpoint`.

# Fields

This type has no public fields.

# Usage

Pass `GIImplicitMidpoint()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIImplicitMidpoint <: GeometricIntegratorAlgorithm end

"""
    GISRK3()

Fixed-step third-order symplectic Runge-Kutta wrapper for `GeometricIntegrators.SRK3`.

# Fields

This type has no public fields.

# Usage

Pass `GISRK3()` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GISRK3 <: GeometricIntegratorAlgorithm end

"""
    GIGLRK(s)

Fixed-step Gauss-Legendre Runge-Kutta wrapper for `GeometricIntegrators.Gauss`.

# Arguments

- `s::Integer`: Number of stages in the Gauss-Legendre method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.Gauss`.

# Usage

Pass `GIGLRK(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIGLRK <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GIRadauIA(s)

Fixed-step Radau IA wrapper for `GeometricIntegrators.RadauIA`.

# Arguments

- `s::Integer`: Number of stages in the Radau IA method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.RadauIA`.

# Usage

Pass `GIRadauIA(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRadauIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GIRadauIIA(s)

Fixed-step Radau IIA wrapper for `GeometricIntegrators.RadauIIA`.

# Arguments

- `s::Integer`: Number of stages in the Radau IIA method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.RadauIIA`.

# Usage

Pass `GIRadauIIA(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GIRadauIIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIA(s)

Fixed-step Lobatto IIIA wrapper for `GeometricIntegrators.LobattoIIIA`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIA method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIA`.

# Usage

Pass `GILobattoIIIA(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIB(s)

Fixed-step Lobatto IIIB wrapper for `GeometricIntegrators.LobattoIIIB`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIB method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIB`.

# Usage

Pass `GILobattoIIIB(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIB <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIC(s)

Fixed-step Lobatto IIIC wrapper for `GeometricIntegrators.LobattoIIIC`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIC method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIC`.

# Usage

Pass `GILobattoIIIC(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIC <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIC̄(s)

Fixed-step Lobatto IIIC-bar wrapper using the Lobatto IIIC tableau from
GeometricIntegrators.jl.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIC-bar method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIC`.

# Usage

Pass `GILobattoIIIC̄(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIC̄ <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIID(s)

Fixed-step Lobatto IIID wrapper for `GeometricIntegrators.LobattoIIID`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIID method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIID`.

# Usage

Pass `GILobattoIIID(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIID <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIE(s)

Fixed-step Lobatto IIIE wrapper for `GeometricIntegrators.LobattoIIIE`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIE method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIE`.

# Usage

Pass `GILobattoIIIE(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIE <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIF(s)

Fixed-step Lobatto IIIF wrapper for `GeometricIntegrators.LobattoIIIF`.

# Arguments

- `s::Integer`: Number of stages in the Lobatto IIIF method.

# Fields

- `s::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIF`.

# Usage

Pass `GILobattoIIIF(s)` to `solve(prob, alg; dt)` as described by
[`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIF <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GISymplecticEulerA()

Fixed-step symplectic Euler A wrapper for `GeometricIntegrators.SymplecticEulerA`.

# Fields

This type has no public fields.

# Usage

Pass `GISymplecticEulerA()` to `solve(prob, alg; dt)` for second-order or dynamical
ODE problems as described by [`GeometricIntegratorAlgorithm`](@ref).
"""
struct GISymplecticEulerA <: GeometricIntegratorAlgorithm end

"""
    GISymplecticEulerB()

Fixed-step symplectic Euler B wrapper for `GeometricIntegrators.SymplecticEulerB`.

# Fields

This type has no public fields.

# Usage

Pass `GISymplecticEulerB()` to `solve(prob, alg; dt)` for second-order or dynamical
ODE problems as described by [`GeometricIntegratorAlgorithm`](@ref).
"""
struct GISymplecticEulerB <: GeometricIntegratorAlgorithm end

"""
    GILobattoIIIAIIIB(n)

Fixed-step Lobatto IIIA-IIIB partitioned wrapper for
`GeometricIntegrators.LobattoIIIAIIIB`.

# Arguments

- `n::Integer`: Number of stages in the partitioned method.

# Fields

- `n::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIAIIIB`.

# Usage

Pass `GILobattoIIIAIIIB(n)` to `solve(prob, alg; dt)` for second-order or dynamical
ODE problems as described by [`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIAIIIB <: GeometricIntegratorAlgorithm
    n::Int
end

"""
    GILobattoIIIBIIIA(n)

Fixed-step Lobatto IIIB-IIIA partitioned wrapper for
`GeometricIntegrators.LobattoIIIBIIIA`.

# Arguments

- `n::Integer`: Number of stages in the partitioned method.

# Fields

- `n::Int`: Number of stages passed to `GeometricIntegrators.LobattoIIIBIIIA`.

# Usage

Pass `GILobattoIIIBIIIA(n)` to `solve(prob, alg; dt)` for second-order or dynamical
ODE problems as described by [`GeometricIntegratorAlgorithm`](@ref).
"""
struct GILobattoIIIBIIIA <: GeometricIntegratorAlgorithm
    n::Int
end
