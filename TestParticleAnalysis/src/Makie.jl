using Makie
using CurrentSheetTestParticle: field_lines, solve_param, field_lines_distance
using CurrentSheetTestParticle: RotationalDiscontinuity

function plot_gc_field_lines!(sol, B; idxs = (1, 2, 3), kwargs...)
    fl0_sol, flf_sol = field_lines(sol, B; kwargs...)
    lines!(fl0_sol; idxs, color = Makie.wong_colors()[3], label = "Initial Field Line")
    lines!(flf_sol; idxs, color = Makie.wong_colors()[4], label = "Final Field Line")
    return fl0_sol, flf_sol
end

function plot_trajectory!(ax, sol; sol_field = missing, sol_gc = missing, add_field_lines = true)

    traj_label(i, μ1, ϕ0) = "T$(i): ϕ₀=$(round(ϕ0, digits = 1))°, cos(α₁)=$(round(μ1, digits = 1))"

    # Plot the actual trajectory
    start_point_coords = sol.prob.u0[1:3] # Get x,y,z of the first point
    scatter!(ax, Point3f(start_point_coords); marker = :star5, markersize = 15, color = Makie.wong_colors()[end - 1]) # Plot the star

    plot!(sol, idxs = (1, 2, 3))

    # Plot the guiding center trajectory if computed
    if sol_gc !== missing
        plot!(sol_gc, idxs = (1, 2, 3), color = :red, label = "Guiding Center")
    end

    # Plot Magnetic Field Line
    if sol_field !== missing
        plot!(sol_field, idxs = (1, 2, 3), label = "Magnetic Field Line", color = :gray, alpha = 0.5)
    end

    if add_field_lines
        idxs = (1, 2, 3)
        fl0_sol, fl1_sol = field_lines(sol; dtmax = 1)
        fl0 = lines!(fl0_sol; idxs, color = :black)
        fl1 = lines!(fl1_sol; idxs, color = :gray)
        axislegend(ax, [fl0, fl1], ["Initial", "Final"], field_line_label(sol), position = (0.9, 0.8))
        # Add the distance between the field lines
        dist = field_lines_asym_distance(sol)
        text!(ax, "dR = $(round(dist, digits = 2))", position = (0.9, 0.7))
    end
    return ax
end

function plot_trajectory(sol; kwargs...)
    f = Figure()
    ax = Axis3(f[1, 1], aspect = :data)
    plot_trajectory!(ax, sol; kwargs...)
    return f
end

field_line_label(B::RotationalDiscontinuity) = L"\text{Field Lines: } θ = %$(B.θ)^∘,\ β = %$(B.β)^∘"
field_line_label(sol) = field_line_label(sol.prob.p.B)
