using QuadGK: quadgk
import OrdinaryDiffEqVerner as DE
using LinearAlgebra: normalize

export trapping_time

"""
Compute the derivative of the position vector along the field line.
drds = B(r) / |B(r)|
"""
field_line_ode(r, B_field, t) = normalize(B_field(r))

function solve_fl(r, B_field, args...; tspan = DEFAULT_TSPAN, kwargs...)
    prob = ODEProblem(field_line_ode, get_r(r), tspan, B_field)
    return solve(prob, Vern9(), args...; verbose = false, kwargs...)
end

# https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#Example-2:-Terminating-an-Integration
function field_lines(gc0, gcf, B; tmax = 100, zmax = nothing, kwargs...)
    zmax = @something zmax maximum(abs ∘ get_z, [gc0, gcf])
    condition(u, _, _) = abs(u[3]) > zmax
    callback = DE.DiscreteCallback(condition, DE.terminate!)

    Bz0 = B(gc0)[3] # this only works for constant Bz
    fl0_sol = solve_fl(gc0, B; tspan = (0, -sign(Bz0 * gc0[3]) * tmax), callback, kwargs...)
    flf_sol = solve_fl(gcf, B; tspan = (0, -sign(Bz0 * gcf[3]) * tmax), callback, kwargs...)
    return fl0_sol, flf_sol
end

function field_lines(sol; kwargs...)
    B = 𝓑(sol)
    gc0 = guiding_center(sol[1], B)
    gcf = guiding_center(sol[end], B)
    zmax = max(maximum(abs ∘ get_z, sol.u), maximum(abs ∘ get_z, [gc0, gcf])) + eps()
    return field_lines(gc0, gcf, B; tmax = 100 * sol.t[end], zmax, kwargs...)
end

"Asymptotic distance between two line solutions."
@views function field_lines_asym_distance(u1s, u2s, B)
    direction = B([0, 0, Inf])
    p1 = argmax(get_z, u1s)
    p2 = argmax(get_z, u2s)
    @assert p1[3] > 5 "p1 is not well above the current sheet (it is not in the asymptotic region)"
    @assert p2[3] > 5 "p2 is not well above the current sheet (it is not in the asymptotic region)"
    distance(p1, p2, direction)
end

function field_lines_min_distance(sol; dtmax = 1.0, kwargs...)
    fl0_sol, flf_sol = field_lines(sol; dtmax, kwargs...)
    return distance(fl0_sol.u, flf_sol.u)
end

# a small `dtmax` is required to accurately compute the minimum distance
function field_lines_distance(sol; kwargs...)
    dR_perp_min = field_lines_min_distance(sol, kwargs...)
    dR_perp_asym = field_lines_asym_distance(sol)
    return (; dR_perp_min, dR_perp_asym)
end

"""
    field_lines_asym_distance(c::RotationalDiscontinuity, r0, r1)

Calculate the shortest distance between the two asymptotic field lines, measured perpendicular to the (asymptotic) magnetic field direction.

The two field lines are defined by the `c::RotationalDiscontinuity` and the two points `r0` and `r1`.
"""
function field_lines_asym_distance(c::RotationalDiscontinuity, r0, r1)
    return norm(field_lines_asym_perp_displacement(c, r0, r1))
end

function field_lines_asym_perp_displacement(c::RotationalDiscontinuity, r0, r1)
    x0, y0, z0 = r0
    x1, y1, z1 = r1
    β = deg2rad(c.β)
    s = c.sign
    tθ = tand(c.θ)

    # evaluate them at the same large z, say z→+∞
    Isin = quadgk(z -> sin(β * tanh(z)), z1, z0)[1]
    Icos = quadgk(z -> cos(β * tanh(z)), z1, z0)[1]
    Δx = (x1 - x0) + tθ * Isin
    Δy = (y1 - y0) + s * tθ * Icos
    Δ𝐫 = SA[Δx, Δy, 0]
    𝐝 = c(Inf)

    Δ𝐫_⊥ = Δ𝐫 - (Δ𝐫 ⋅ 𝐝) * 𝐝 / (𝐝 ⋅ 𝐝)
    return Δ𝐫_⊥
end


function field_lines_asym_distance(sol)
    B = 𝓑(sol)
    gc0 = guiding_center(sol[1], B)
    gc1 = guiding_center(sol[end], B)
    return field_lines_asym_distance(B, gc0, gc1)
end

function field_lines_asym_distance(sn::SolutionNorm)
    return field_lines_asym_distance(sn.sol) * L₀(sn)
end

"""
    projection(𝐱, sol)

Reference: https://www.wikiwand.com/en/articles/Projection_(linear_algebra)

A simple case occurs when the orthogonal projection is onto a line. If ``𝐮`` is a unit vector on the line, then the projection is given by the outer product ``𝐮𝐮^T``.
"""
function projection(𝐱, sol)
    B = 𝓑(sol)
    gc1 = guiding_center(sol[end], B)
    return 𝐱 + field_lines_asym_perp_displacement(B, 𝐱, gc1)
end

"""
    approach_time(𝐮, 𝐁)

Calculate the time for the guiding center of the particle to reach current sheet (z=0) accounting for 
field variation along trajectory, assuming no scattering.
"""
function approach_time(𝐮, 𝐁; z0 = 0.0, rtol = 1.0e-4)
    𝐫 = guiding_center(𝐮, 𝐁)
    zi = get_z(𝐫)
    v_para = v_parallel(𝐮, 𝐁)
    # Integrate time along field line from zi to 0
    # dt = ds / v_∥ where ds is field line element
    # For field line: dr/ds = B̂
    # So dz/ds = B_z / |B|
    # Therefore: dt = dz / (v_∥ * B_z / |B|)
    integrand(z) = 1 / normalize(𝐁(z))[3]

    # Integrate from current z to 0
    s, err = quadgk(integrand, zi, z0; rtol)
    return s / v_para
end

leaving_time(𝐮, 𝐁; kw...) = -approach_time(𝐮, 𝐁; kw...)

# A very simple model of trapping time
function trapping_time(sol; kw...)
    𝐁 = 𝓑(sol)
    t0 = approach_time(sol.u[1], 𝐁; kw...)
    tf = leaving_time(sol.u[end], 𝐁; kw...)
    return sol.t[end] - t0 - tf
end

trapping_time(sn::SolutionNorm; kw...) = trapping_time(sn.sol; kw...) * T₀(sn)
