function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractODEProblem{uType, tType, isinplace},
        alg::AlgType,
        timeseries = nothing, ts = nothing, ks = nothing;
        verbose = true,
        save_start = true, dt = nothing,
        timeseries_errors = true,
        callback = nothing, alias_u0 = false,
        kwargs...
    ) where {
        uType, tType, isinplace,
        AlgType <: GeometricIntegratorAlgorithm,
    }
    if dt === nothing
        error("dt required for fixed timestep methods.")
    end

    isstiff = !(
        alg isa Union{
            GIImplicitEuler, GIImplicitMidpoint,
            GISRK3, GIGLRK, GIRadauIA, GIRadauIIA,
        }
    )

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(prob.f isa DiffEqBase.AbstractParameterizedFunction) && isstiff
            if DiffEqBase.has_tgrad(prob.f)
                @warn "Explicit t-gradient given to this stiff solver is ignored."
                warned = true
            end
            if DiffEqBase.has_jac(prob.f)
                @warn "Explicit Jacobian given to this stiff solver is ignored."
                warned = true
            end
        end
        warned && warn_compat()
    end

    if callback !== nothing
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

    _alg = get_method_from_alg(alg)
    needs_solver = requires_newton_solver(alg)

    if prob.problem_type isa DiffEqBase.StandardODEProblem
        # Create function wrapper for GeometricIntegrators API
        # GeometricIntegrators expects: v(v, t, q, params)
        # DiffEqBase provides: f(du, u, p, t) for inplace or f(u, p, t) for out-of-place
        if !isinplace && u isa AbstractArray
            v! = (v, t, q, params) -> (v .= vec(prob.f(reshape(q, sizeu), p, t)); nothing)
        elseif !(u isa Vector{Float64})
            v! = (
                v, t, q,
                params,
            ) -> (prob.f(reshape(v, sizeu), reshape(q, sizeu), p, t); nothing)
        else
            v! = (v, t, q, params) -> prob.f(v, q, p, t)
        end

        ode = GeometricIntegrators.ODEProblem(v!, prob.tspan, dt, vec(prob.u0))
        if needs_solver
            sol = integrate(ode, _alg; solver = Newton())
        else
            sol = integrate(ode, _alg)
        end

        if save_start
            start_idx = 1
            ts = prob.tspan[1]:dt:prob.tspan[end]
        else
            start_idx = 2
            ts = (prob.tspan[1] + dt):dt:prob.tspan[end]
        end

        if u0 isa DiffEqBase.RecursiveArrayTools.ArrayPartition
            _timeseries = [copy(vec(q)) for q in sol.q]
        elseif u0 isa AbstractArray
            _timeseries = [reshape(copy(vec(q)), sizeu) for q in sol.q]
        else
            _timeseries = [q[1] for q in sol.q]
        end

    elseif prob.problem_type isa DiffEqBase.AbstractDynamicalODEProblem
        # For second order / partitioned problems
        # PODE expects: v(v, t, q, p, params), f(f, t, q, p, params)
        # DiffEqBase SecondOrderODEProblem has f.f1 and f.f2:
        # For inplace: f.f1(dv, v, u, p, t) - for acceleration (dv/dt = f(v, u))
        # For out-of-place: f.f1(v, u, p, t) returns dv

        v! = (v, t, q, p_state, params) -> (v .= p_state)  # dq/dt = p

        # Handle both inplace and out-of-place problems
        if isinplace
            f! = (
                f_out, t, q, p_state,
                params,
            ) -> (prob.f.f1.f(f_out, p_state, q, p, t); nothing)  # dp/dt = f1(p, q)
        else
            f! = (
                f_out, t, q, p_state,
                params,
            ) -> (f_out .= prob.f.f1.f(p_state, q, p, t); nothing)
        end

        pode = GeometricIntegrators.PODEProblem(
            v!, f!, prob.tspan, dt, vec(prob.u0.x[1]), vec(prob.u0.x[2])
        )
        if needs_solver
            sol = integrate(pode, _alg; solver = Newton())
        else
            sol = integrate(pode, _alg)
        end

        if save_start
            start_idx = 1
            ts = prob.tspan[1]:dt:prob.tspan[end]
        else
            start_idx = 2
            ts = (prob.tspan[1] + dt):dt:prob.tspan[end]
        end

        _timeseries = [
            DiffEqBase.RecursiveArrayTools.ArrayPartition(copy(vec(q)), copy(vec(p)))
                for (q, p) in zip(sol.q, sol.p)
        ]
    end

    return DiffEqBase.build_solution(
        prob, alg, ts, _timeseries[start_idx:end],
        timeseries_errors = timeseries_errors,
        retcode = ReturnCode.Success
    )
end

# Use multiple dispatch for requires_newton_solver - matches get_method_from_alg pattern
# Algorithms that require a Newton solver for implicit stages
requires_newton_solver(::GIRadauIA) = true
requires_newton_solver(::GIRadauIIA) = true
requires_newton_solver(::GILobattoIIIA) = true
requires_newton_solver(::GILobattoIIIC) = true
requires_newton_solver(::GILobattoIIIC̄) = true
requires_newton_solver(::GILobattoIIID) = true
requires_newton_solver(::GILobattoIIIE) = true
requires_newton_solver(::GILobattoIIIF) = true
# Default: no Newton solver required
requires_newton_solver(::GeometricIntegratorAlgorithm) = false

# Use multiple dispatch for better type inference and cleaner code
get_method_from_alg(::GIEuler) = ExplicitEuler()
get_method_from_alg(::GIMidpoint) = ExplicitMidpoint()
get_method_from_alg(::GIHeun2) = Heun2()
get_method_from_alg(::GIHeun3) = Heun3()
get_method_from_alg(::GIRalston2) = Ralston2()
get_method_from_alg(::GIRalston3) = Ralston3()
get_method_from_alg(::GIRunge) = Runge2()
get_method_from_alg(::GIKutta) = Kutta3()
get_method_from_alg(::GIRK4) = RK4()
get_method_from_alg(::GIRK416) = RK416()
get_method_from_alg(::GIRK438) = RK438()
get_method_from_alg(::GISSPRK3) = SSPRK3()
get_method_from_alg(::GICrankNicolson) = CrankNicolson()
get_method_from_alg(::GIKraaijevangerSpijker) = KraaijevangerSpijker()
get_method_from_alg(::GIQinZhang) = QinZhang()
get_method_from_alg(::GICrouzeix) = Crouzeix()
get_method_from_alg(::GIImplicitEuler) = ImplicitEuler()
get_method_from_alg(::GIImplicitMidpoint) = ImplicitMidpoint()
get_method_from_alg(::GISRK3) = SRK3()
get_method_from_alg(alg::GIGLRK) = Gauss(alg.s)
get_method_from_alg(alg::GIRadauIA) = RadauIA(alg.s)
get_method_from_alg(alg::GIRadauIIA) = RadauIIA(alg.s)
get_method_from_alg(alg::GILobattoIIIA) = LobattoIIIA(alg.s)
get_method_from_alg(alg::GILobattoIIIB) = LobattoIIIB(alg.s)
get_method_from_alg(alg::GILobattoIIIC) = LobattoIIIC(alg.s)
get_method_from_alg(alg::GILobattoIIIC̄) = LobattoIIIC(alg.s)  # LobattoIIIC̄ not available as standalone
get_method_from_alg(alg::GILobattoIIID) = LobattoIIID(alg.s)
get_method_from_alg(alg::GILobattoIIIE) = LobattoIIIE(alg.s)
get_method_from_alg(alg::GILobattoIIIF) = LobattoIIIF(alg.s)
get_method_from_alg(::GISymplecticEulerA) = SymplecticEulerA()
get_method_from_alg(::GISymplecticEulerB) = SymplecticEulerB()
get_method_from_alg(alg::GILobattoIIIAIIIB) = LobattoIIIAIIIB(alg.n)
get_method_from_alg(alg::GILobattoIIIBIIIA) = LobattoIIIBIIIA(alg.n)
