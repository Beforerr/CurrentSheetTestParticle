module CurrentSheetTestParticle
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using SciMLBase: isinplace
using LinearAlgebra
using StaticArrays
using UnPack
using Moshi.Match: @match

export RotationalDiscontinuity, TD_B_field
export solve_params
export w_ϕ_pairs, init_state, init_states, init_states_pm, filter_wϕs!
export isoutofdomain_params
export ProblemParamsBase, ProblemParams, RDProblemParams
export trace_normalized_B!, trace_normalized_B

include("utils.jl")
include("field.jl")
include("types.jl")
include("state.jl")
include("fieldline.jl")
include("pa.jl")
include("equations.jl")

abstol = 1e-7 # Defaults to 1e-5
reltol = 1e-7 # Defaults to 1e-3
maxiters = 1e6 # Defaults to 1e5
dtmin = 1e-4 # Reduce the computation for domain checking (as `isoutofdomain` will reject and reducd time step until a step is accepted)
const DEFAULT_SOLVER = AutoVern9(Rodas4P())
const DEFAULT_DIFFEQ_KWARGS = (; abstol, reltol, maxiters, dtmin)
const DEFAULT_BORIS_KWARGS = (; dt=1e-2, savestepinterval=1)
const DEFAULT_TSPAN = (0, 256)
const ez = SA[0, 0, 1]

@kwdef struct ProblemParamsBase{F,V,A,I,T,D}
    B::F
    v::V = 1
    alg::A = :AutoVern9
    init_kwargs::I = (; Nw=8, Nϕ=8)
    tspan::T = DEFAULT_TSPAN
    diffeq::D = DEFAULT_DIFFEQ_KWARGS
end

function RDProblemParams(; θ=DEFAULT_θ, β=DEFAULT_β, sign=DEFAULT_SIGN, kwargs...)
    Bf = RotationalDiscontinuity(; θ, β, sign)
    ProblemParamsBase(; B=Bf, kwargs...)
end

const ProblemParams = RDProblemParams

isoutofdomain_z(z_max) = (u, p, t) -> abs(u[3]) > abs(z_max)

function isoutofdomain_params(v)
    z_init = init_z_pos(v)
    z_max = 2 * z_init
    return isoutofdomain_z(z_max)
end

"""
Solve the system of ODEs.
"""
function solve_params(B, u0s::AbstractVector, tspan=DEFAULT_TSPAN, f=trace_normalized_B, iipv::Val{iip}=Val(false); E=E0func, alg=DEFAULT_SOLVER, diffeq=DEFAULT_DIFFEQ_KWARGS, kwargs...) where {iip}
    param = Parameters(; E, B)
    u0s = iip ? u0s : [SVector{6}(u0) for u0 in u0s]
    ensemble_solve(f, u0s, tspan, param, _alg(alg), iipv; merge(diffeq, kwargs)...)
end

"""
Type-stable version of `solve` for `solve_params`.
"""
function ensemble_solve(f, u0s::AbstractVector, tspan, param, alg, ::Val{iip}; kwargs...) where {iip}
    # Making tspan float for Type stability: https://discourse.julialang.org/t/type-instability-when-solving-odeproblem/105687
    prob = ODEProblem{iip}(f, u0s[1], Float64.(tspan), param)
    prob_func = (prob, i, repeat=nothing) -> remake(prob, u0=u0s[i])
    ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
    solve(ensemble_prob, alg, EnsembleThreads(); trajectories=length(u0s), kwargs...)
end

function solve_params(d::ProblemParamsBase; kwargs...)
    @unpack B, v, alg, init_kwargs, diffeq = d
    u0s, wϕs = init_states_pm(B, v; init_kwargs...)

    isoutofdomain = isoutofdomain_params(v)
    sol = solve_params(B, u0s, d.tspan; alg, diffeq, isoutofdomain, kwargs...)
    return sol, (wϕs, B)
end
end
