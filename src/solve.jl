function RDProblemParams(; θ = DEFAULT_θ, β = DEFAULT_β, sign = DEFAULT_SIGN, kwargs...)
    Bf = RotationalDiscontinuity(; θ, β, sign)
    return ProblemParamsBase(; B = Bf, kwargs...)
end

struct OutOfDomainCondition
    z_max::Float64
end

(o::OutOfDomainCondition)(u, _, _) = abs(u[3]) > o.z_max

"""
Solve the system of ODEs.
"""
function solve_params(B, u0s, tspan = DEFAULT_TSPAN, f = trace_normalized_B; E = F0func, alg = DEFAULT_SOLVER, diffeq = DEFAULT_DIFFEQ_KWARGS, kwargs...)
    param = Parameters(; E, B)
    return ensemble_solve(f, u0s, tspan, param, _alg(alg), Val(false); merge(diffeq, kwargs)...)
end

function solve_param(B, u0, tspan = DEFAULT_TSPAN, f = trace_normalized_B, iipv::Val{iip} = Val(false); alg = DEFAULT_SOLVER, diffeq = DEFAULT_DIFFEQ_KWARGS, kwargs...) where {iip}
    u0 = iip ? u0 : SVector{6, Float64}(u0)
    prob = ODEProblem{false}(f, u0, Float64.(tspan), B)
    return solve(prob, _alg(alg); merge(diffeq, kwargs)...)
end

"""
Type-stable version of `solve` for `solve_params`.
"""
function ensemble_solve(f, u0s::AbstractVector, tspan, param, alg, ::Val{iip}; ensemblealg = EnsembleThreads(), kwargs...) where {iip}
    # Making tspan float for Type stability: https://discourse.julialang.org/t/type-instability-when-solving-odeproblem/105687
    prob = ODEProblem{iip}(f, u0s[1], Float64.(tspan), param)
    prob_func = (prob, i, repeat = nothing) -> remake(prob, u0 = u0s[i])
    ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)
    return solve(ensemble_prob, alg, ensemblealg; trajectories = length(u0s), kwargs...)
end

function solve_params(d::ProblemParamsBase; load = false, kwargs...)
    @unpack B, v, alg, init_kwargs, diffeq = d
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    condition = OutOfDomainCondition(abs(2 * init_z_pos(v)))
    callback = DE.DiscreteCallback(condition, DE.terminate!; save_positions = (false, false))

    path = joinpath(".cache", string(hash(d), ".jld2"))
    sol = produce_or_load(path; load) do
        solve_params(B, u0s, d.tspan, trace_normalized_B; alg, diffeq, callback, kwargs...)
    end
    return sol, (wϕs, B, u0s)
end

𝓑(sol::ODESolution) = 𝓑(sol.prob.p)
𝓑(s::SolutionNorm) = 𝓑(s.sol)
