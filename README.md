# GeometricIntegratorsDiffEq.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/GeometricIntegratorsDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/GeometricIntegratorsDiffEq.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/GeometricIntegratorsDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/GeometricIntegratorsDiffEq.jl?branch=master)

This package contains bindings for GeometricIntegrators.jl to allow it to be used with the
JuliaDiffEq common interface. For more information on using the solvers from this
package, see the [DifferentialEquations.jl documentation](https://docs.sciml.ai/stable/).

## Common API Usage

This library adds the common interface to GeometricIntegrators.jl's solvers. [See the DifferentialEquations.jl documentation for details on the interface](https://docs.sciml.ai/stable/index.html). Following the Lorenz example from [the ODE tutorial](https://docs.sciml.ai/stable/tutorials/ode_example/), we can solve this using `GIEuler` via the following:

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

The options available in `solve` are documented [at the common solver options page](https://docs.sciml.ai/stable/basics/common_solver_opts/). The available methods are documented [at the ODE solvers page](https://docs.sciml.ai/stable/solvers/ode_solve#GeometricIntegrators.jl-1).
