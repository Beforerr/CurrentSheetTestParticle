module Analysis
using ..CurrentSheetTestParticle: cos_pitch_angle
using DataFrames, DataFramesMeta

export process_sols, process!

function extract_info(sol)
    u0 = sol.prob.u0
    u1 = sol.u[end]
    t1 = sol.t[end]
    return (; u0, u1, t1)
end

process_sols(sols) = DataFrame(map(extract_info, sols))

function process_sols(sols, B, wϕs)
    df = process_sols(sols)
    df.wϕ0 .= wϕs
    df.B .= B
    return df
end

function process!(df)
    "wϕ0" in names(df) && @chain df begin
        @rtransform!(:μ0 = :wϕ0[1], :ϕ0 = :wϕ0[2], :μ1 = cos_pitch_angle(:u1, :B))
        @transform!(:μ0 = clamp.(:μ0, -1, 1), :μ1 = clamp.(:μ1, -1, 1))
        @transform!(:α0 = acosd.(:μ0), :α1 = acosd.(:μ1))
        @transform!(:Δα = :α1 .- :α0, :Δμ = :μ1 .- :μ0)
    end

    return df
end

end