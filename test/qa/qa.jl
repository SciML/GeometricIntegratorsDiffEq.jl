using SciMLTesting, GeometricIntegratorsDiffEq, Test

run_qa(
    GeometricIntegratorsDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # The handful of names below are still non-public in their owners' latest
        # registered releases (SciMLBase 3.24.0, DiffEqBase 7.5.7): a few SciMLBase
        # API-extension points / problem types that have not yet been `public`-declared,
        # and `ArrayPartition` reached through the non-public `DiffEqBase.RecursiveArrayTools`
        # module. Ignore them until they go public upstream.
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractParameterizedFunction,  # SciMLBase
                :StandardODEProblem,             # SciMLBase
                :__solve,                        # SciMLBase
                :has_tgrad,                      # SciMLBase
                :unwrapped_f,                    # SciMLBase
                :RecursiveArrayTools,            # DiffEqBase (module path to ArrayPartition)
            ),
        ),
    ),
)
