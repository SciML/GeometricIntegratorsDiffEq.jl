using SciMLTesting, GeometricIntegratorsDiffEq, Test

run_qa(
    GeometricIntegratorsDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # The names below are still non-public in their owner's latest registered
        # release (verified on Julia 1.12 against SciMLBase 3.28.1 that
        # `Base.ispublic` returns `false` for each): SciMLBase API-extension
        # points / problem types not yet `public`-declared. Ignore them until they
        # go public upstream.
        all_qualified_accesses_are_public = (;
            ignore = (
                :StandardODEProblem,             # SciMLBase
                :__solve,                        # SciMLBase
                :unwrapped_f,                    # SciMLBase
            ),
        ),
    ),
)
