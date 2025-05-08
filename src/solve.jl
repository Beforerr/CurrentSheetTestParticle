function RDProblemParams(; θ=DEFAULT_θ, β=DEFAULT_β, sign=DEFAULT_SIGN, kwargs...)
    Bf = RotationalDiscontinuity(; θ, β, sign)
    ProblemParamsBase(; B=Bf, kwargs...)
end

function isoutofdomain_params(v)
    z_init = init_z_pos(v)
    z_max = 2 * z_init
    return OutOfDomainZ(abs(z_max))
end

struct OutOfDomainZ
    z_max::Float64
end

(o::OutOfDomainZ)(u, p, t) = abs(u[3]) > o.z_max

"""
Solve the system of ODEs.
"""
function solve_params(B, u0s, tspan=DEFAULT_TSPAN, f=trace_normalized_B, iipv::Val{iip}=Val(false); E=F0func, alg=DEFAULT_SOLVER, diffeq=DEFAULT_DIFFEQ_KWARGS, kwargs...) where {iip}
    param = Parameters(; E, B)
    u0s = iip ? u0s : [SVector{6}(u0) for u0 in u0s]
    ensemble_solve(f, u0s, tspan, param, _alg(alg), iipv; merge(diffeq, kwargs)...)
end

"""
Type-stable version of `solve` for `solve_params`.
"""
function ensemble_solve(f, u0s::AbstractVector, tspan, param, alg, ::Val{iip}; ensemblealg=EnsembleThreads(), kwargs...) where {iip}
    # Making tspan float for Type stability: https://discourse.julialang.org/t/type-instability-when-solving-odeproblem/105687
    prob = ODEProblem{iip}(f, u0s[1], Float64.(tspan), param)
    prob_func = (prob, i, repeat=nothing) -> remake(prob, u0=u0s[i])
    ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
    solve(ensemble_prob, alg, ensemblealg; trajectories=length(u0s), kwargs...)
end

function solve_params(d::ProblemParamsBase; kwargs...)
    @unpack B, v, alg, init_kwargs, diffeq = d
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    isoutofdomain = isoutofdomain_params(v)
    sol = solve_params(B, u0s, d.tspan; alg, diffeq, isoutofdomain, kwargs...)
    return sol, (wϕs, B)
end
