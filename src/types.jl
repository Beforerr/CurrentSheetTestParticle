const abstol = 1.0e-7 # Defaults to 1e-5
const reltol = 1.0e-7 # Defaults to 1e-3
const maxiters = 1.0e6 # Defaults to 1e5
const DEFAULT_SOLVER = AutoVern9(Rodas4P())
const DEFAULT_TSPAN = (0.0, 256.0)
const DEFAULT_DIFFEQ_KWARGS = (; abstol, reltol, maxiters)

F0func(::Any) = SVector(0.0, 0.0, 0.0)

@kwdef struct Parameters{TE, TB, TF, Q, M}
    E::TE = F0func
    B::TB = F0func
    F::TF = F0func
    q::Q = 1.0
    m::M = 1.0
end

@kwdef struct ProblemParamsBase{F, V, A, I, T, D}
    B::F
    v::V = 1
    alg::A = :AutoVern9
    init_kwargs::I = (; Nw = 8, Nϕ = 8)
    tspan::T = DEFAULT_TSPAN
    diffeq::D = DEFAULT_DIFFEQ_KWARGS
end

struct SolutionNorm{S, N}
    sol::S
    norm::N
end

@inline function Base.getproperty(sol::SolutionNorm, name::Symbol)
    name in (:sol, :norm) && return getfield(sol, name)
    return getproperty(sol.sol, name)
end

L₀(sn::SolutionNorm) = sn.norm.L₀
T₀(sn::SolutionNorm) = sn.norm.T₀
V₀(sn::SolutionNorm) = sn.norm.V₀
