# GeometricIntegratorsDiffEq.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/GeometricIntegratorsDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/GeometricIntegratorsDiffEq.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/GeometricIntegratorsDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl?branch=master)

This package contains bindings for GeometricIntegrators.jl to allow it to be used with the
JuliaDiffEq common interface. For more information on using the solvers from this
package, see the [DifferentialEquations.jl documentation](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/).

## Installation

Since GeometricIntegrators.jl is not registered, it must be loaded separately. Note that GeometricIntegrators.jl 
segfaults on non-Linux machines, and thus GeometricIntegratorsDiffEq.jl will not work on non-Linux as well.

## Common API Usage

This library adds the common interface to GeometricIntegrators.jl's solvers. [See the DifferentialEquations.jl documentation for details on the interface](http://docs.juliadiffeq.org/latest/index.html). Following the Lorenz example from [the ODE tutorial](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html), we can solve this using `GIEuler` via the following:

```julia
using GeometricIntegratorsDiffEq
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob,GIEuler(),dt=0.1)
using Plots; plot(sol,vars=(1,2,3))
```

The options available in `solve` are documented [at the common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html). The available methods are documented [at the ODE solvers page](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html#GeometricIntegrators.jl-1).
