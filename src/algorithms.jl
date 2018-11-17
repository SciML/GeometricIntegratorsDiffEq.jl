abstract type GeometricIntegratorAlgorithm <: DiffEqBase.DEAlgorithm end
struct GIEuler <: GeometricIntegratorAlgorithm end
struct GIMidpoint <: GeometricIntegratorAlgorithm end
struct GIHeun <: GeometricIntegratorAlgorithm end
struct GIKutta <: GeometricIntegratorAlgorithm end
struct GIERK4 <: GeometricIntegratorAlgorithm end
struct GIERK438 <: GeometricIntegratorAlgorithm end
struct GIImplicitEuler <: GeometricIntegratorAlgorithm end
struct GIImplicitMidpoint <: GeometricIntegratorAlgorithm end
struct GIRadIIA2 <: GeometricIntegratorAlgorithm end
struct GIRadIIA3 <: GeometricIntegratorAlgorithm end
struct GISRK3 <: GeometricIntegratorAlgorithm end
struct GIGLRK <: GeometricIntegratorAlgorithm
    s::Int
end
struct GISymplecticEulerA <: GeometricIntegratorAlgorithm end
struct GISymplecticEulerB <: GeometricIntegratorAlgorithm end
struct GILobattoIIIAIIIB2 <: GeometricIntegratorAlgorithm end
struct GILobattoIIIBIIIA2 <: GeometricIntegratorAlgorithm end
