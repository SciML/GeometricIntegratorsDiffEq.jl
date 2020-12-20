abstract type GeometricIntegratorAlgorithm <: DiffEqBase.DEAlgorithm end
struct GIEuler <: GeometricIntegratorAlgorithm end
struct GIMidpoint <: GeometricIntegratorAlgorithm end
struct GIHeun2 <: GeometricIntegratorAlgorithm end
struct GIHeun3 <: GeometricIntegratorAlgorithm end
struct GIRalston2 <: GeometricIntegratorAlgorithm end
struct GIRalston3 <: GeometricIntegratorAlgorithm end
struct GIRunge <: GeometricIntegratorAlgorithm end
struct GIKutta <: GeometricIntegratorAlgorithm end
struct GIRK4 <: GeometricIntegratorAlgorithm end
struct GIRK416 <: GeometricIntegratorAlgorithm end
struct GIRK438 <: GeometricIntegratorAlgorithm end
struct GISSPRK3 <: GeometricIntegratorAlgorithm end
struct GICrankNicolson <: GeometricIntegratorAlgorithm end
struct GIKraaijevangerSpijker <: GeometricIntegratorAlgorithm end
struct GIQinZhang <: GeometricIntegratorAlgorithm end
struct GICrouzeix <: GeometricIntegratorAlgorithm end
struct GIImplicitEuler <: GeometricIntegratorAlgorithm end
struct GIImplicitMidpoint <: GeometricIntegratorAlgorithm end
struct GISRK3 <: GeometricIntegratorAlgorithm end
struct GIGLRK <: GeometricIntegratorAlgorithm
    s::Int
end
struct GIRadauIA <: GeometricIntegratorAlgorithm
    s::Int
end
struct GIRadauIIA <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIA <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIB <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIC <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIC̄ <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIID <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIE <: GeometricIntegratorAlgorithm
    s::Int
end
struct GILobattoIIIF <: GeometricIntegratorAlgorithm
    s::Int
end
struct GISymplecticEulerA <: GeometricIntegratorAlgorithm end
struct GISymplecticEulerB <: GeometricIntegratorAlgorithm end
struct GILobattoIIIAIIIB2 <: GeometricIntegratorAlgorithm end
struct GILobattoIIIBIIIA2 <: GeometricIntegratorAlgorithm end
