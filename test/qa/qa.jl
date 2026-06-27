using SciMLTesting, GeometricIntegratorsDiffEq, Test

run_qa(
    GeometricIntegratorsDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # The names below are still non-public in their owners' latest registered
        # releases (verified on Julia 1.12 against SciMLBase 3.27.0, DiffEqBase 7.6.0:
        # `Base.ispublic` returns `false` for each): a few SciMLBase API-extension
        # points / problem types not yet `public`-declared, and `ArrayPartition`
        # reached through the non-public `DiffEqBase.RecursiveArrayTools` module.
        # Ignore them until they go public upstream.
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
