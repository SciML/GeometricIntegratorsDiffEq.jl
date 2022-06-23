function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType, tType, isinplace},
                            alg::AlgType,
                            timeseries = [], ts = [], ks = [];
                            verbose = true,
                            save_start = true, dt = nothing,
                            timeseries_errors = true,
                            callback = nothing, alias_u0 = false,
                            kwargs...) where {uType, tType, isinplace,
                                              AlgType <: GeometricIntegratorAlgorithm}
    if dt == nothing
        error("dt required for fixed timestep methods.")
    end

    N = ceil(Int, (prob.tspan[end] - prob.tspan[1]) / dt)

    isstiff = !(typeof(alg) <: Union{GIImplicitEuler, GIImplicitMidpoint,
                      GISRK3, GIGLRK, GIRadauIA, GIRadauIIA})

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

    if !isinplace && typeof(u) <: AbstractArray
        f! = (t, u, du) -> (du[:] = vec(prob.f(reshape(u, sizeu), p, t)); nothing)
    elseif !(typeof(u) <: Vector{Float64})
        f! = (t, u, du) -> (prob.f(reshape(du, sizeu), reshape(u, sizeu), p, t); nothing)
    elseif typeof(prob.problem_type) <: DiffEqBase.StandardODEProblem
        f! = (t, u, du) -> prob.f(du, u, p, t)
    end

    _alg = get_tableau_from_alg(alg)
    if typeof(prob.problem_type) <: DiffEqBase.StandardODEProblem
        ode = ODE(f!, vec(prob.u0))
    elseif typeof(prob.problem_type) <: DiffEqBase.AbstractDynamicalODEProblem
        ode = PODE((t, u, v, du) -> prob.f.f2(du, v, u, t),
                   (t, u, v, dv) -> prob.f.f1(dv, v, u, t),
                   vec(prob.u0.x[1]), vec(prob.u0.x[2]))
    end
    sol = integrate(ode, _alg, dt, N)

    if save_start
        start_idx = 1
        ts = prob.tspan[1]:dt:prob.tspan[end]
    else
        start_idx = 2
        ts = (prob.tspan[1] + dt):dt:prob.tspan[end]
    end

    if typeof(u0) <: DiffEqBase.RecursiveArrayTools.ArrayPartition
        _timeseries = sol.q.d
    elseif typeof(u0) <: Union{AbstractArray}
        _timeseries = map(x -> reshape(x, sizeu), sol.q.d)
    else
        _timeseries = vec(sol.q)
    end

    DiffEqBase.build_solution(prob, alg, ts, _timeseries,
                              timeseries_errors = timeseries_errors,
                              retcode = :Success)
end

function get_tableau_from_alg(alg)
    typeof(alg) == GIEuler && (_alg = TableauExplicitEuler())
    typeof(alg) == GIMidpoint && (_alg = TableauExplicitMidpoint())
    typeof(alg) == GIHeun2 && (_alg = TableauHeun2())
    typeof(alg) == GIHeun3 && (_alg = TableauHeun3())
    typeof(alg) == GIRalston2 && (_alg = TableauRalston2())
    typeof(alg) == GIRalston3 && (_alg = TableauRalston3())
    typeof(alg) == GIRunge && (_alg = TableauRunge())
    typeof(alg) == GIKutta && (_alg = TableauKutta())
    typeof(alg) == GIRK4 && (_alg = TableauRK4())
    typeof(alg) == GIRK416 && (_alg = TableauRK416())
    typeof(alg) == GIRK438 && (_alg = TableauRK438())
    typeof(alg) == GISSPRK3 && (_alg = TableauSSPRK3())
    typeof(alg) == GICrankNicolson && (_alg = TableauCrankNicolson())
    typeof(alg) == GIKraaijevangerSpijker && (_alg = TableauKraaijevangerSpijker())
    typeof(alg) == GIQinZhang && (_alg = TableauQinZhang())
    typeof(alg) == GICrouzeix && (_alg = TableauCrouzeix())
    typeof(alg) == GIImplicitEuler && (_alg = TableauImplicitEuler())
    typeof(alg) == GIImplicitMidpoint && (_alg = TableauImplicitMidpoint())
    typeof(alg) == GISRK3 && (_alg = TableauSRK3())
    typeof(alg) == GIGLRK && (_alg = TableauGauss(alg.s))
    typeof(alg) == GIRadauIA && (_alg = TableauRadauIA(alg.s))
    typeof(alg) == GIRadauIIA && (_alg = TableauRadauIIA(alg.s))
    typeof(alg) == GILobattoIIIA && (_alg = TableauLobattoIIIA(alg.s))
    typeof(alg) == GILobattoIIIB && (_alg = TableauLobattoIIIB(alg.s))
    typeof(alg) == GILobattoIIIC && (_alg = TableauLobattoIIIC(alg.s))
    typeof(alg) == GILobattoIIIC̄ && (_alg = TableauLobattoIIIC̄(alg.s))
    typeof(alg) == GILobattoIIID && (_alg = TableauLobattoIIID(alg.s))
    typeof(alg) == GILobattoIIIE && (_alg = TableauLobattoIIIE(alg.s))
    typeof(alg) == GILobattoIIIF && (_alg = TableauLobattoIIIF(alg.s))
    typeof(alg) == GISymplecticEulerA && (_alg = TableauLobattoIIIAIIIB(2))
    typeof(alg) == GISymplecticEulerB && (_alg = TableauLobattoIIIBIIIA(2))
    typeof(alg) == GILobattoIIIAIIIB && (_alg = TableauLobattoIIIAIIIB(2))
    typeof(alg) == GILobattoIIIBIIIA && (_alg = TableauLobattoIIIBIIIA(2))
    _alg
end
