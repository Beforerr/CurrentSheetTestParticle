"""
Compute the derivative of the position vector along the field line.
drds = B(r) / |B(r)|
"""
field_line_ode(r, B_field, t) = normalize(B_field(r))

function solve_fl(r, B_field, args...; tspan=DEFAULT_TSPAN, kwargs...)
    prob = ODEProblem(field_line_ode, r, tspan, B_field)
    return solve(prob, Vern9(), args...; verbose=false, kwargs...)
end

function field_lines(sol, B; kwargs...)
    gc0 = guiding_center(sol[1], B)
    gcf = guiding_center(sol[end], B)
    zmax = maximum(abs ∘ get_z, sol.u)
    isoutofdomain = (u, p, t) -> abs(get_z(u)) > zmax
    Bz0 = B(gc0)[3] # this only works for constant Bz
    tmax = 100 * sol.t[end]
    fl0_sol = solve_fl(gc0, B; tspan=(0, -sign(Bz0 * gc0[3]) * tmax), isoutofdomain, kwargs...)
    flf_sol = solve_fl(gcf, B; tspan=(0, -sign(Bz0 * gcf[3]) * tmax), isoutofdomain, kwargs...)
    return fl0_sol, flf_sol
end

"""
Asymptotic distance between two line solutions
"""
@views function field_lines_asym_distance(u1s, u2s, B)
    direction = B([0, 0, Inf])
    p1 = argmax(get_z, u1s)
    p2 = argmax(get_z, u2s)
    distance(p1, p2, direction)
end

"""
Field lines distances
"""
function field_lines_distance(sol1, sol2, B)
    dR_perp_min = distance(get_r.(sol1.u), get_r.(sol2.u))
    dR_perp_asym = field_lines_asym_distance(sol1.u, sol2.u, B)
    return (; dR_perp_min, dR_perp_asym)
end

for sym in [:field_lines_distance, :field_lines_asym_distance]
    @eval $sym(sol, B) = $sym(field_lines(sol, B)..., B)
end