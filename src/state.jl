const ez = SA[0, 0, 1]
const DEFAULT_Z_INIT_0 = 5
# TODO: check the initial position effect
"""
Particle with `pos=1/-1` is initialized above/below the current sheet
"""
init_z_pos(v; pos=-1, z_init_0=DEFAULT_Z_INIT_0) = pos * (abs(z_init_0) + 2 * abs(v))
init_pos(v; kw...) = [0, 0, init_z_pos(v; kw...)]

function init_state(v, w, ϕ, B::Function; kw...)
    r₀ = init_pos(v; kw...)
    v₀ = init_v(v, w, ϕ, B(r₀))
    return vcat(r₀, v₀)
end

function init_states(B::Function, v, wϕs::Vector; kw...)
    return map(wϕs) do wϕ
        init_state(v, wϕ..., B; kw...)
    end
end

function init_states_pm(B::Function, v; kw...)
    wϕs = w_ϕ_pairs(; kw...)

    wϕs_below = filter(wϕ -> wϕ[1] > 0, wϕs)
    u0s_below = init_states(B, v, wϕs_below; pos=-1)

    wϕs_above = filter(wϕ -> wϕ[1] < 0, wϕs)
    u0s_above = init_states(B, v, wϕs_above; pos=1)
    return vcat(u0s_below, u0s_above), vcat(wϕs_below, wϕs_above)
end

function local_B_coord(B; e1=ez)
    e_para = B |> normalize  # parallel to B
    e_perp1 = e1 × e_para |> normalize
    e_perp2 = e_para × e_perp1 |> normalize
    return e_para, e_perp1, e_perp2
end

"""
    init_v(v, w, ϕ, ...)

Initialize the velocity vector of the particle with magnitude `v`, cosine pitch angle `w`, and azimuthal angle `ϕ` at position `r`.
"""
function init_v(v, w, ϕ)
    v_para = v * w
    v_perp1 = v * sqrt(1 - w^2) * cos(ϕ)
    v_perp2 = v * sqrt(1 - w^2) * sin(ϕ)
    return [v_para, v_perp1, v_perp2]
end

function init_v(v, w, ϕ, e_para, e_perp1, e_perp2)
    return [e_para e_perp1 e_perp2] * init_v(v, w, ϕ)
end

init_v(v, w, ϕ, B0; e1=ez) = init_v(v, w, ϕ, local_B_coord(B0; e1=e1)...)
init_v(v, w, ϕ, r, B; e1=ez) = init_v(v, w, ϕ, B(r); e1=e1)

"""
    inverse_v(v1, v2, v3, ...)

Inverse the velocity vector (v1, v2, v3) to get the magnitude `v`, cosine pitch angle `w`, and azimuthal angle `ϕ`.
Note: we assume ϕ is in the range of [0, 2π).
"""
function inverse_v(v1, v2, v3)
    v = sqrt(v1^2 + v2^2 + v3^2)
    w = v1 / v
    ϕ = rem2pi(atan(v3, v2), RoundDown)
    return v, w, ϕ
end

function inverse_v(v1, v2, v3, e_para, e_perp1, e_perp2)
    v_local = [e_para e_perp1 e_perp2] \ [v1, v2, v3]
    return inverse_v(v_local...)
end

inverse_v(v, B0; e1=ez) = inverse_v(v..., local_B_coord(B0; e1)...)
inverse_v(u, B::Function; kw...) = @views inverse_v(u[4:6], B(u[1:3]); kw...)

init_range(::Nothing; start, stop, length=2) = range(; start, stop, length)
init_range(x::AbstractArray; kw...) = x
init_range(x::Number; kw...) = [x]

"""
Generate a grid of w and phi pairs for the particles.
"""
function w_ϕ_pairs(; Nw=8, Nϕ=8, w=nothing, ϕ=nothing)
    ws = init_range(w; start=-1, stop=1, length=Nw)
    ϕs = init_range(ϕ; start=0, stop=2π - 2π / Nϕ, length=Nϕ)
    w_ϕ_pairs(ws, ϕs)
end

w_ϕ_pairs(ws::AbstractArray, ϕs::AbstractArray) = [(w, ϕ) for w in ws for ϕ in ϕs]

"""
Filter pitch angle corresponding to particles moving toward the current sheet
"""
function filter_wϕs!(wϕs, θ)
    if θ == pi / 2
        return wϕs
    elseif θ < pi / 2
        return filter!(wϕ -> wϕ[1] > 0, wϕs)
    elseif θ > pi / 2
        return filter!(wϕ -> wϕ[1] < 0, wϕs)
    end
end