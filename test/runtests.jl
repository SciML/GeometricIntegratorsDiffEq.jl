using GeometricIntegratorsDiffEq
using Base.Test

using DiffEqProblemLibrary
using GeometricIntegrators

prob = prob_ode_2Dlinear

sol = solve(prob,GIEuler(),dt=0.1)
sol = solve(prob,GIMidpoint(),dt=0.1)
sol = solve(prob,GIHeun(),dt=0.1)
sol = solve(prob,GIKutta(),dt=0.1)
sol = solve(prob,GIERK4(),dt=0.1)
sol = solve(prob,GIERK438(),dt=0.1)
sol = solve(prob,GIImplicitEuler(),dt=0.1)
sol = solve(prob,GIImplicitMidpoint(),dt=0.1)
sol = solve(prob,GIRadIIA2(),dt=0.1)
sol = solve(prob,GIRadIIA3(),dt=0.1)
sol = solve(prob,GISRK3(),dt=0.1)
sol = solve(prob,GIGLRK(2),dt=0.1)

u0 = zeros(2)
v0 = ones(2)
f2 = function (t,u,v,dv)
  dv .= -u
end
function (::typeof(f2))(::Type{Val{:analytic}}, x, y0)
  u0, v0 = y0
  ArrayPartition(u0*cos(x) + v0*sin(x), -u0*sin(x) + v0*cos(x))
end
prob = SecondOrderODEProblem(f2,u0,v0,(0.0,5.0))

sol = solve(prob,GISymplecticEulerA(),dt=0.1)
sol = solve(prob,GISymplecticEulerB(),dt=0.1)
sol = solve(prob,GILobattoIIIAIIIB2(),dt=0.1)
sol = solve(prob,GILobattoIIIBIIIA2(),dt=0.1)
