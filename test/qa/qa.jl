using SciMLTesting, GeometricIntegratorsDiffEq, Test

run_qa(
    GeometricIntegratorsDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # All qualified accesses are sourced from each name's owner, but a handful of
        # SciMLBase API-extension points / problem types are non-public (not exported
        # or `public`-declared) and `ArrayPartition` is reached through the non-public
        # `DiffEqBase.RecursiveArrayTools` module. These go public as the base libraries
        # release; ignore them until then. Source package noted per name.
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractDynamicalODEProblem,    # SciMLBase
                :AbstractODEAlgorithm,           # SciMLBase
                :AbstractODEProblem,             # SciMLBase
                :AbstractParameterizedFunction,  # SciMLBase
                :StandardODEProblem,             # SciMLBase
                :__solve,                        # SciMLBase
                :build_solution,                 # SciMLBase
                :has_jac,                        # SciMLBase
                :has_tgrad,                      # SciMLBase
                :unwrapped_f,                    # SciMLBase
                :Success,                        # SciMLBase.ReturnCode (enum member; flagged only on Julia 1.10)
                :RecursiveArrayTools,            # DiffEqBase (module path to ArrayPartition)
            ),
        ),
    ),
)
