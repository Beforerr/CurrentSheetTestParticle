const DEFAULT_θ = 45
const DEFAULT_β = 90
const DEFAULT_SIGN = 1

abstract type AbstractField <: Function end
abstract type MagneticField <: AbstractField end

struct ConstField{T} <: MagneticField
    𝐁::T
end

(c::ConstField)(args...) = c.𝐁

"""
    RotationalDiscontinuity

Rotating magnetic field. 

# Arguments
- `θ` : azimuthal angle, angle between Bn and B0
- `β` : half of rotation angle
- `sign` : sign of rotation, 1 for left-handed, -1 for right-handed

# Notes
φ = β * tanh(z) is the polar angle
"""
@kwdef struct RotationalDiscontinuity{T1,T2,T3} <: MagneticField
    B::T1 = 1
    θ::T2 = DEFAULT_θ
    β::T3 = DEFAULT_β
    sign::Int = DEFAULT_SIGN
end

"""
# Keyword Arguments
- `dir` : direction where the field depends on, 1, 2, or 3
"""
B(r, conf::RotationalDiscontinuity; dir=3) = conf(r; dir)

function (c::RotationalDiscontinuity)(z)
    B = c.B
    θ = c.θ
    φ = c.β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    By = c.sign * B * sind(θ) * cosd(φ)
    return SVector(Bx, By, Bz)
end

(c::RotationalDiscontinuity)(r::AbstractArray; dir=3) = c(r[dir])
(c::RotationalDiscontinuity)(r, t; dir=3) = c(r; dir)


function TD_B_field(r; dir=3, B=1, By=0, θ=DEFAULT_θ, β=DEFAULT_β, kw...)
    z = r[dir]
    φ = β * tanh(z)
    Bz = B * cosd(θ)
    Bx = B * sind(θ) * sind(φ)
    return SVector(Bx, By, Bz)
end

TD_B_field(; kw...) = r -> TD_B_field(r; kw...)