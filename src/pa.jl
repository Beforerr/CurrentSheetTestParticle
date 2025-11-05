using LinearAlgebra: norm, dot

export cos_pitch_angle, cos_pitch_angles, pa
export magnetic_moment, magnetic_moments

"Cos pitch angle"
function cos_pitch_angle(u, B)
    v = u[SA[4, 5, 6]]
    return clamp(v ⋅ B / (norm(v) * norm(B)), -1, 1)
end

cos_pitch_angle(u, B::Function) = cos_pitch_angle(u, B(u))

function cos_pitch_angles(sol)
    return @. cos_pitch_angle(sol.u, 𝓑(sol))
end

const pa = cos_pitch_angle
