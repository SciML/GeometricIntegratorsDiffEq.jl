using Documenter, GeometricIntegratorsDiffEq

makedocs(
    modules = [GeometricIntegratorsDiffEq],
    sitename = "GeometricIntegratorsDiffEq.jl",
    pages = ["Home" => "index.md"],
    checkdocs = :exports,
    doctest = true,
)
