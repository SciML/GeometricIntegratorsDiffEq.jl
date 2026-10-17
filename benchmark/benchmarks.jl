using GeometricIntegratorsDiffEq, BenchmarkTools

const SUITE = BenchmarkGroup()

# Linear oscillator
function f_linear!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
    return nothing
end
prob = ODEProblem(f_linear!, [1.0, 0.0], (0.0, 50.0))

# Dynamical (second-order) problem: q'' = -q
function dyn_f1!(dv, v, q, p, t)
    dv[1] = -q[1]
    return nothing
end
function dyn_f2!(dq, v, q, p, t)
    dq[1] = v[1]
    return nothing
end
dyn_prob = DynamicalODEProblem(
    DynamicalODEFunction(dyn_f1!, dyn_f2!), [0.0], [1.0], (0.0, 50.0)
)

# =============================================================================
# Explicit geometric integrators
# =============================================================================

SUITE["explicit"] = BenchmarkGroup()

SUITE["explicit"]["GIEuler"] = @benchmarkable solve($prob, GIEuler(); dt = 0.1)
SUITE["explicit"]["GIMidpoint"] = @benchmarkable solve($prob, GIMidpoint(); dt = 0.1)
SUITE["explicit"]["GIHeun3"] = @benchmarkable solve($prob, GIHeun3(); dt = 0.1)
SUITE["explicit"]["GIRK4"] = @benchmarkable solve($prob, GIRK4(); dt = 0.1)
SUITE["explicit"]["GIRK438"] = @benchmarkable solve($prob, GIRK438(); dt = 0.1)

# =============================================================================
# Implicit geometric integrators
# =============================================================================

SUITE["implicit"] = BenchmarkGroup()

SUITE["implicit"]["GIImplicitMidpoint"] = @benchmarkable solve(
    $prob, GIImplicitMidpoint(); dt = 0.1
)
SUITE["implicit"]["GIGLRK2"] = @benchmarkable solve($prob, GIGLRK(2); dt = 0.1)

# =============================================================================
# Dynamical problem (partitioned structure)
# =============================================================================

SUITE["dynamical"] = BenchmarkGroup()

SUITE["dynamical"]["GISymplecticEulerA"] = @benchmarkable solve(
    $dyn_prob, GISymplecticEulerA(); dt = 0.1
)
