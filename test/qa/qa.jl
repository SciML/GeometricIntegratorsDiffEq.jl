using SciMLTesting, GeometricIntegratorsDiffEq, Test

# The SciML common interface GeometricIntegratorsDiffEq deliberately reexports so that
# `using GeometricIntegratorsDiffEq` is enough to build and solve a problem, as the
# README and the docstring examples do. Owned and documented upstream; kept in sync with
# the reexport `export` block in src/GeometricIntegratorsDiffEq.jl.
const REEXPORTS = (
    :DynamicalODEFunction, :DynamicalODEProblem, :EnsembleAnalysis, :EnsembleDistributed,
    :EnsembleProblem, :EnsembleSerial, :EnsembleSolution, :EnsembleSplitThreads,
    :EnsembleSummary, :EnsembleThreads, :NullParameters, :ODEFunction, :ODEProblem,
    :ODESolution, :ReturnCode, :SecondOrderODEProblem, :remake, :solve,
    :successful_retcode,
)

run_qa(GeometricIntegratorsDiffEq; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using
    # GeometricIntegratorsDiffEq`, so the allow-list cannot drift into approving names
    # the package no longer provides.
    @testset "$name" for name in REEXPORTS
        @test name in names(GeometricIntegratorsDiffEq)
        @test isdefined(@__MODULE__, name)
    end
end
