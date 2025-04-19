const abstol = 1e-7 # Defaults to 1e-5
const reltol = 1e-7 # Defaults to 1e-3
const maxiters = 1e6 # Defaults to 1e5
const dtmin = 1e-4 # Reduce the computation for domain checking (as `isoutofdomain` will reject and reducd time step until a step is accepted)
const DEFAULT_SOLVER = AutoVern9(Rodas4P())
const DEFAULT_TSPAN = (0., 256.)
const DEFAULT_DIFFEQ_KWARGS = (; abstol, reltol, maxiters, dtmin)

F0func(x) = SVector(0.0, 0.0, 0.0)

@kwdef struct Parameters{TE,TB,TF,Q,M}
    E::TE = F0func
    B::TB = F0func
    F::TF = F0func
    q::Q = 1.0
    m::M = 1.0
end

Bfunc(t::Tuple) = t[3]
Bfunc(p::Parameters) = p.B


@kwdef struct ProblemParamsBase{F,V,A,I,T,D}
    B::F
    v::V = 1
    alg::A = :AutoVern9
    init_kwargs::I = (; Nw=8, Nϕ=8)
    tspan::T = DEFAULT_TSPAN
    diffeq::D = DEFAULT_DIFFEQ_KWARGS
end
