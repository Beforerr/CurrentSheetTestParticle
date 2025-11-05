module Analysis
using ..CurrentSheetTestParticle: cos_pitch_angle, magnetic_moment, field_lines_asym_distance, 𝓑, get_r, get_v
using LinearAlgebra: norm
using DataFrames, DataFramesMeta

export extract_info, process_sols, process!

function extract_info(sol)
    u0 = sol.prob.u0
    u1 = sol.u[end]
    t1 = sol.t[end]
    tmax = sol.prob.tspan[2]
    # dR_perp_asym = field_lines_asym_distance(sol)
    return (; u0, u1, t1, tmax)
end

function process_sols(sols)
    df = DataFrame(map(extract_info, sols))
    df.B .= map(𝓑, sols)
    return @rtransform! df @astable begin
        :μ0 = cos_pitch_angle(:u0, :B)
        :μ1 = cos_pitch_angle(:u1, :B)
        :α0 = acosd.(:μ0)
        :α1 = acosd.(:μ1)
        :M0 = magnetic_moment(:u0, :B)
        :M1 = magnetic_moment(:u1, :B)
        :M_max = norm(get_v(:u0))^2 / (2 * norm((:B)(get_r(:u0))))
        :ΔM_norm = (:M1 - :M0) / :M_max
    end
end

function process_sols(prob, sols)
    df = process_sols(sols)
    df.v .= prob.v
    return process!(df)
end

function process!(df)
    "dR_perp_asym" in names(df) && @transform!(df, :dR_perp_asym_norm = :dR_perp_asym ./ :v)
    @transform!(df, :leave = :t1 .!= :tmax)
    return df
end

end
