module CurrentSheetTestParticle
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner: AutoVern9, @unpack
using LinearAlgebra
using StaticArrays
using Moshi.Match: @match

export RotationalDiscontinuity, TD_B_field
export solve_params
export w_ϕ_pairs, init_state, init_states, init_states_pm, filter_wϕs!
export isoutofdomain_params
export ProblemParamsBase, RDProblemParams
export trace_normalized_B!, trace_normalized_B

include("utils.jl")
include("field.jl")
include("types.jl")
include("state.jl")
include("fieldline.jl")
include("pa.jl")
include("equations.jl")
include("solve.jl")
end
