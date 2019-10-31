function DiffEqBase.__solve(
    prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
    alg::AlgType,
    timeseries=[],ts=[],ks=[];
    verbose=true,
    save_start=true, dt = nothing,
    timeseries_errors=true,
    callback=nothing, alias_u0=false, kwargs...) where {uType,tType,isinplace,AlgType<:GeometricIntegratorAlgorithm}

    if dt == nothing
        error("dt required for fixed timestep methods.")
    end

    N = ceil(Int,(prob.tspan[end]-prob.tspan[1]) / dt)

    isstiff = !(typeof(alg) <: Union{GIImplicitEuler,GIImplicitMidpoint,GIRadIIA2,
                                     GIRadIIA3,GISRK3,GIGLRK})

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(typeof(prob.f) <: DiffEqBase.AbstractParameterizedFunction) && isstiff
            if DiffEqBase.has_tgrad(prob.f)
                warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
            if DiffEqBase.has_jac(prob.f)
                warn("Explicit Jacobian given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    if callback != nothing
        error("GeometricIntegrators is not compatible with callbacks.")
    end

    u0 = prob.u0

    if typeof(u0) <: Number
        u = [u0]
    else
        if alias_u0
            u = u0
        else
            u = deepcopy(u0)
        end
    end

    if typeof(u) <: Tuple
        sizeu = size(u[1])
    else
        sizeu = size(u)
    end

    p = prob.p

    if !isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = vec(prob.f(reshape(u,sizeu),p,t)); nothing)
    elseif !(typeof(u)<:Vector{Float64})
        f! = (t,u,du) -> (prob.f(reshape(du,sizeu),reshape(u,sizeu),p,t); nothing)
    elseif typeof(prob.problem_type) <: DiffEqBase.StandardODEProblem
        f! = (t,u,du) -> prob.f(du,u,p,t)
    end

    _alg = get_tableau_from_alg(alg)
    if typeof(prob.problem_type) <: DiffEqBase.StandardODEProblem
        ode = ODE(f!, vec(prob.u0))
    elseif typeof(prob.problem_type) <: DiffEqBase.AbstractDynamicalODEProblem
        ode = PODE((t,u,v,du)->prob.f.f2(du,v,u,t),
                   (t,u,v,dv)->prob.f.f1(dv,v,u,t),
                   vec(prob.u0.x[1]),vec(prob.u0.x[2]))
    end
    integrator = Integrator(ode,_alg,dt)
    sol = integrate(integrator, N)

    if save_start
        start_idx = 1
        ts = prob.tspan[1]:dt:prob.tspan[end]
    else
        start_idx = 2
        ts = (prob.tspan[1]+dt):dt:prob.tspan[end]
    end

    if typeof(u0) <: DiffEqBase.RecursiveArrayTools.ArrayPartition
        _timeseries = Vector{typeof(u0.x[1])}(undef,0)
        for i=start_idx:sol.q.nt
            push!(_timeseries, reshape(view(sol.q, :, i-1)', length(u0.x[1])))
        end
    elseif typeof(u0) <: Union{AbstractArray}
        _timeseries = Vector{uType}(undef,0)
        for i=start_idx:sol.q.nt
            push!(_timeseries, reshape(view(sol.q, :, i-1)', sizeu))
        end
    else
        _timeseries = vec(sol.q)
    end

    DiffEqBase.build_solution(prob,  alg, ts, _timeseries,
                   timeseries_errors = timeseries_errors,
                   retcode = :Success)
end

function get_tableau_from_alg(alg)
    typeof(alg) == GIEuler && (_alg = getTableauExplicitEuler())
    typeof(alg) == GIMidpoint && (_alg = getTableauExplicitMidpoint())
    typeof(alg) == GIHeun && (_alg = getTableauHeun())
    typeof(alg) == GIKutta && (_alg = getTableauKutta())
    typeof(alg) == GIERK4 && (_alg = getTableauERK4())
    typeof(alg) == GIERK438 && (_alg = getTableauERK438())
    typeof(alg) == GIImplicitEuler && (_alg = getTableauImplicitEuler())
    typeof(alg) == GIImplicitMidpoint && (_alg = getTableauImplicitMidpoint())
    typeof(alg) == GIRadIIA2 && (_alg = getTableauRadIIA2())
    typeof(alg) == GIRadIIA3 && (_alg = getTableauRadIIA3())
    typeof(alg) == GISRK3 && (_alg = getTableauSRK3())
    typeof(alg) == GIGLRK && (_alg = getTableauGLRK(alg.s))
    typeof(alg) == GISymplecticEulerA && (_alg = getTableauSymplecticEulerA())
    typeof(alg) == GISymplecticEulerB && (_alg = getTableauSymplecticEulerB())
    typeof(alg) == GILobattoIIIAIIIB2 && (_alg = getTableauLobattoIIIAIIIB2())
    typeof(alg) == GILobattoIIIBIIIA2 && (_alg = getTableauLobattoIIIBIIIA2())
    _alg
end
