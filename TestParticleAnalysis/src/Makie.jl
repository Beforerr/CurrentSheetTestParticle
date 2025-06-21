using Makie
using CurrentSheetTestParticle: field_lines

function plot_gc_field_lines!(sol, B; idxs=(1, 2, 3), kwargs...)
    fl0_sol, flf_sol = field_lines(sol, B; kwargs...)
    lines!(fl0_sol; idxs, color=Makie.wong_colors()[3], label="Initial Field Line")
    lines!(flf_sol; idxs, color=Makie.wong_colors()[4], label="Final Field Line")
    return fl0_sol, flf_sol
end