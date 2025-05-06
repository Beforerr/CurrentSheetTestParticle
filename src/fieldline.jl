"""
Compute the derivative of the position vector along the field line.
drds = B(r) / |B(r)|
"""
function field_line_ode!(drds, r, p, t)
    B_field = p[1]
    drds .= normalize(B_field(r))
end

function solve_fl(r, B_field, args...; tspan=DEFAULT_TSPAN, kwargs...)
    prob = ODEProblem(field_line_ode!, r, tspan, (B_field,))
    return solve(prob, Vern9(), args...; verbose=false, kwargs...)
end

function field_lines(sol, B; kwargs...)
    gc0 = guiding_center(sol[1], B)
    gcf = guiding_center(sol[end], B)
    isoutofdomain = (u, p, t) -> abs(u[3]) > maximum(abs.(sol[3, :]))
    Bz0 = B(gc0)[3] # this only works for constant Bz
    tmax = 100 * sol.t[end]
    fl0_sol = solve_fl(gc0, B; tspan=(0, -sign(Bz0 * gc0[3]) * tmax), isoutofdomain, kwargs...)
    flf_sol = solve_fl(gcf, B; tspan=(0, -sign(Bz0 * gcf[3]) * tmax), isoutofdomain, kwargs...)
    return fl0_sol, flf_sol
end

"""
Distance of two line solutions
"""
@views distance(sol1::ODESolution, sol2::ODESolution) = distance(sol1[1:3, :], sol2[1:3, :])

"""
Asymptotic distance between two line solutions
"""

@views function field_lines_asym_distance(u1s, u2s, B)
    direction = B([0, 0, Inf])
    p1 = u1s[argmax(u1s[3, :])]
    p2 = u2s[argmax(u2s[3, :])]
    distance(p1, p2, direction)
end

field_lines_asym_distance(sol1::ODESolution, sol2::ODESolution, B) =
    field_lines_asym_distance(sol1.u, sol2.u, B)

"""
Field lines distances
"""
function field_lines_distance(sol1, sol2, B)
    dR_perp_min = distance(sol1, sol2)
    dR_perp_asym = field_lines_asym_distance(sol1, sol2, B)
    return (; dR_perp_min, dR_perp_asym)
end

for sym in [:field_lines_distance, :field_lines_asym_distance]
    @eval $sym(sol, B) = $sym(field_lines(sol, B)..., B)
end

# Accelerated version of `field_lines_asym_distance`
function _field_lines_asym_distance(sol, B)
    direction = B([0, 0, Inf])
    gc0 = guiding_center(sol[1], B)
    gcf = guiding_center(sol[end], B)
    if sign(gc0[3]) == sign(gcf[3])
        return distance(gc0, gcf, direction)
    else
        isoutofdomain = (u, p, t) -> abs(u[3]) > maximum(abs.(sol[3, :]))
        Bz0 = B(gc0)[3] # this only works for constant Bz
        tmax = 100 * sol.t[end]
        flf_sol = solve_fl(gcf, B; tspan=(0, -sign(Bz0 * gcf[3]) * tmax), isoutofdomain, save_everystep=false)
        return distance(gc0, flf_sol[end], direction)
    end
end