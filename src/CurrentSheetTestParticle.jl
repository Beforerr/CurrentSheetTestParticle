module CurrentSheetTestParticle
using JLD2
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner: AutoVern9, Vern9
using LinearAlgebra
using StaticArrays
using Moshi.Match: @match
using UnPack: @unpack
using Unitful: q, mp, @u_str, upreferred
using Unitful: Velocity, Mass, BField, Density, Charge, Energy
using Magnetostatics: MagneticField, RotationalDiscontinuity

export 𝓑, MagneticField, RotationalDiscontinuity
export solve_param, solve_params
export w_ϕ_pairs, init_state, init_states, init_states_pm, filter_wϕs!
export ProblemParamsBase, RDProblemParams
export trace_normalized_B!, trace_normalized_B
export guiding_center, guiding_centers, field_lines_asym_distance
export field_lines
export energy, velocity

const DEFAULT_Z_INIT_0 = 5

# 𝓑 is the mapping of parameters to fields 𝓑 : p ↦ 𝐁(𝐫)
function 𝓑 end

include("field.jl")
include("types.jl")
include("utils.jl")
include("state.jl")
include("fieldline.jl")
include("pa.jl")
include("equations.jl")
include("solve.jl")
include("analysis.jl")

is_outside(u; z_init_0 = DEFAULT_Z_INIT_0) = abs(u[3]) > z_init_0


end
