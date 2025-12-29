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

    isstiff = !(alg isa Union{GIImplicitEuler, GIImplicitMidpoint,
        GISRK3, GIGLRK, GIRadauIA, GIRadauIIA})

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(prob.f isa DiffEqBase.AbstractParameterizedFunction) && isstiff
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

    if u0 isa Number
        u = [u0]
    else
        if alias_u0
            u = u0
        else
            u = deepcopy(u0)
        end
    end

    if u isa Tuple
        sizeu = size(u[1])
    else
        sizeu = size(u)
    end

    p = prob.p

    if !isinplace && u isa AbstractArray
        f! = (t, u, du) -> (du[:] = vec(prob.f(reshape(u, sizeu), p, t)); nothing)
    elseif !(u isa Vector{Float64})
        f! = (t, u, du) -> (prob.f(reshape(du, sizeu), reshape(u, sizeu), p, t); nothing)
    elseif prob.problem_type isa DiffEqBase.StandardODEProblem
        f! = (t, u, du) -> prob.f(du, u, p, t)
    end

    _alg = get_tableau_from_alg(alg)
    if prob.problem_type isa DiffEqBase.StandardODEProblem
        ode = ODE(f!, vec(prob.u0))
    elseif prob.problem_type isa DiffEqBase.AbstractDynamicalODEProblem
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

    if u0 isa DiffEqBase.RecursiveArrayTools.ArrayPartition
        _timeseries = sol.q.d
    elseif u0 isa AbstractArray
        # sol.q.d is a matrix where rows are state variables and columns are time points
        # Extract each column, skip first column if not saving start
        _timeseries = [reshape(view(sol.q.d, :, i), sizeu) for i in start_idx:size(sol.q.d, 2)]
    else
        _timeseries = vec(sol.q)
    end

    DiffEqBase.build_solution(prob, alg, ts, _timeseries,
        timeseries_errors = timeseries_errors,
        retcode = :Success)
end

# Use dispatch instead of sequential typeof checks for better performance
get_tableau_from_alg(::GIEuler) = TableauExplicitEuler()
get_tableau_from_alg(::GIMidpoint) = TableauExplicitMidpoint()
get_tableau_from_alg(::GIHeun2) = TableauHeun2()
get_tableau_from_alg(::GIHeun3) = TableauHeun3()
get_tableau_from_alg(::GIRalston2) = TableauRalston2()
get_tableau_from_alg(::GIRalston3) = TableauRalston3()
get_tableau_from_alg(::GIRunge) = TableauRunge()
get_tableau_from_alg(::GIKutta) = TableauKutta()
get_tableau_from_alg(::GIRK4) = TableauRK4()
get_tableau_from_alg(::GIRK416) = TableauRK416()
get_tableau_from_alg(::GIRK438) = TableauRK438()
get_tableau_from_alg(::GISSPRK3) = TableauSSPRK3()
get_tableau_from_alg(::GICrankNicolson) = TableauCrankNicolson()
get_tableau_from_alg(::GIKraaijevangerSpijker) = TableauKraaijevangerSpijker()
get_tableau_from_alg(::GIQinZhang) = TableauQinZhang()
get_tableau_from_alg(::GICrouzeix) = TableauCrouzeix()
get_tableau_from_alg(::GIImplicitEuler) = TableauImplicitEuler()
get_tableau_from_alg(::GIImplicitMidpoint) = TableauImplicitMidpoint()
get_tableau_from_alg(::GISRK3) = TableauSRK3()
get_tableau_from_alg(alg::GIGLRK) = TableauGLRK(alg.s)
get_tableau_from_alg(alg::GIRadauIA) = TableauRadauIA(alg.s)
get_tableau_from_alg(alg::GIRadauIIA) = TableauRadauIIA(alg.s)
get_tableau_from_alg(alg::GILobattoIIIA) = TableauLobattoIIIA(alg.s)
get_tableau_from_alg(alg::GILobattoIIIB) = TableauLobattoIIIB(alg.s)
get_tableau_from_alg(alg::GILobattoIIIC) = TableauLobattoIIIC(alg.s)
get_tableau_from_alg(alg::GILobattoIIIC̄) = TableauLobattoIIIC̄(alg.s)
get_tableau_from_alg(alg::GILobattoIIID) = TableauLobattoIIID(alg.s)
get_tableau_from_alg(alg::GILobattoIIIE) = TableauLobattoIIIE(alg.s)
get_tableau_from_alg(alg::GILobattoIIIF) = TableauLobattoIIIF(alg.s)
get_tableau_from_alg(::GISymplecticEulerA) = TableauSymplecticEulerA()
get_tableau_from_alg(::GISymplecticEulerB) = TableauSymplecticEulerB()
get_tableau_from_alg(::GILobattoIIIAIIIB) = TableauLobattoIIIAIIIB2()
get_tableau_from_alg(::GILobattoIIIBIIIA) = TableauLobattoIIIBIIIA2()
