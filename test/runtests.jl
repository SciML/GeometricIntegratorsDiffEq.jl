using SafeTestsets
using Test
using SciMLTesting

run_tests(;
    core = function ()
        @safetestset "GeometricIntegratorsDiffEq" begin
            include("core_tests.jl")
        end
        return @safetestset "Explicit Imports" begin
            include("explicit_imports_test.jl")
        end
    end,
    groups = Dict(
        # NoPre runs the JET and AllocCheck tests in their own environment. The
        # original dispatcher ran these only for GROUP=="NoPre" (never as part of
        # "All"), so NoPre is an env-bearing group kept out of the curated `all`.
        "NoPre" => (;
            env = joinpath(@__DIR__, "NoPre"),
            body = function ()
                @time @safetestset "JET Tests" begin
                    include(joinpath(@__DIR__, "NoPre", "jet.jl"))
                end
                return @time @safetestset "AllocCheck Tests" begin
                    include(joinpath(@__DIR__, "NoPre", "alloc_tests.jl"))
                end
            end,
        ),
    ),
    # The original runtests.jl ran the Core body for GROUP=All and GROUP=Core, and
    # ran NoPre only for GROUP=NoPre (never under "All"). Curate "All" to Core only.
    all = ["Core"],
)
