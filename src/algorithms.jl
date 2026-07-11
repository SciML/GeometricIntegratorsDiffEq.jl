"""
    GeometricIntegratorAlgorithm

Abstract base type for GeometricIntegratorsDiffEq ODE algorithms.
"""
abstract type GeometricIntegratorAlgorithm <: SciMLBase.AbstractODEAlgorithm end

"""
    GIEuler()

GeometricIntegratorsDiffEq wrapper for the explicit Euler method.
"""
struct GIEuler <: GeometricIntegratorAlgorithm end

"""
    GIMidpoint()

GeometricIntegratorsDiffEq wrapper for the explicit midpoint method.
"""
struct GIMidpoint <: GeometricIntegratorAlgorithm end

"""
    GIHeun2()

GeometricIntegratorsDiffEq wrapper for Heun's second-order method.
"""
struct GIHeun2 <: GeometricIntegratorAlgorithm end

"""
    GIHeun3()

GeometricIntegratorsDiffEq wrapper for Heun's third-order method.
"""
struct GIHeun3 <: GeometricIntegratorAlgorithm end

"""
    GIRalston2()

GeometricIntegratorsDiffEq wrapper for Ralston's second-order method.
"""
struct GIRalston2 <: GeometricIntegratorAlgorithm end

"""
    GIRalston3()

GeometricIntegratorsDiffEq wrapper for Ralston's third-order method.
"""
struct GIRalston3 <: GeometricIntegratorAlgorithm end

"""
    GIRunge()

GeometricIntegratorsDiffEq wrapper for Runge's second-order method.
"""
struct GIRunge <: GeometricIntegratorAlgorithm end

"""
    GIKutta()

GeometricIntegratorsDiffEq wrapper for Kutta's third-order method.
"""
struct GIKutta <: GeometricIntegratorAlgorithm end

"""
    GIRK4()

GeometricIntegratorsDiffEq wrapper for the classical fourth-order Runge-Kutta method.
"""
struct GIRK4 <: GeometricIntegratorAlgorithm end

"""
    GIRK416()

GeometricIntegratorsDiffEq wrapper for the RK416 method.
"""
struct GIRK416 <: GeometricIntegratorAlgorithm end

"""
    GIRK438()

GeometricIntegratorsDiffEq wrapper for the RK438 method.
"""
struct GIRK438 <: GeometricIntegratorAlgorithm end

"""
    GISSPRK3()

GeometricIntegratorsDiffEq wrapper for the third-order SSPRK method.
"""
struct GISSPRK3 <: GeometricIntegratorAlgorithm end

"""
    GICrankNicolson()

GeometricIntegratorsDiffEq wrapper for the Crank-Nicolson method.
"""
struct GICrankNicolson <: GeometricIntegratorAlgorithm end

"""
    GIKraaijevangerSpijker()

GeometricIntegratorsDiffEq wrapper for the Kraaijevanger-Spijker method.
"""
struct GIKraaijevangerSpijker <: GeometricIntegratorAlgorithm end

"""
    GIQinZhang()

GeometricIntegratorsDiffEq wrapper for the Qin-Zhang method.
"""
struct GIQinZhang <: GeometricIntegratorAlgorithm end

"""
    GICrouzeix()

GeometricIntegratorsDiffEq wrapper for the Crouzeix method.
"""
struct GICrouzeix <: GeometricIntegratorAlgorithm end

"""
    GIImplicitEuler()

GeometricIntegratorsDiffEq wrapper for the implicit Euler method.
"""
struct GIImplicitEuler <: GeometricIntegratorAlgorithm end

"""
    GIImplicitMidpoint()

GeometricIntegratorsDiffEq wrapper for the implicit midpoint method.
"""
struct GIImplicitMidpoint <: GeometricIntegratorAlgorithm end

"""
    GISRK3()

GeometricIntegratorsDiffEq wrapper for the third-order symplectic Runge-Kutta method.
"""
struct GISRK3 <: GeometricIntegratorAlgorithm end

"""
    GIGLRK(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Gauss-Legendre Runge-Kutta method.
"""
struct GIGLRK <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GIRadauIA(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Radau IA method.
"""
struct GIRadauIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GIRadauIIA(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Radau IIA method.
"""
struct GIRadauIIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIA(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIA method.
"""
struct GILobattoIIIA <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIB(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIB method.
"""
struct GILobattoIIIB <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIC(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIC method.
"""
struct GILobattoIIIC <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIC̄(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIC-bar method.
"""
struct GILobattoIIIC̄ <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIID(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIID method.
"""
struct GILobattoIIID <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIE(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIE method.
"""
struct GILobattoIIIE <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GILobattoIIIF(s)

GeometricIntegratorsDiffEq wrapper for an `s`-stage Lobatto IIIF method.
"""
struct GILobattoIIIF <: GeometricIntegratorAlgorithm
    s::Int
end

"""
    GISymplecticEulerA()

GeometricIntegratorsDiffEq wrapper for the symplectic Euler A method.
"""
struct GISymplecticEulerA <: GeometricIntegratorAlgorithm end

"""
    GISymplecticEulerB()

GeometricIntegratorsDiffEq wrapper for the symplectic Euler B method.
"""
struct GISymplecticEulerB <: GeometricIntegratorAlgorithm end

"""
    GILobattoIIIAIIIB(n)

GeometricIntegratorsDiffEq wrapper for a Lobatto IIIA-IIIB partitioned method.
"""
struct GILobattoIIIAIIIB <: GeometricIntegratorAlgorithm
    n::Int
end

"""
    GILobattoIIIBIIIA(n)

GeometricIntegratorsDiffEq wrapper for a Lobatto IIIB-IIIA partitioned method.
"""
struct GILobattoIIIBIIIA <: GeometricIntegratorAlgorithm
    n::Int
end
