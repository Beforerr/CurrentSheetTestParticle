using CurrentSheetTestParticle: RotationalDiscontinuity
using SpecialFunctions: sinint, cosint
using QuadGK: quadgk

export fieldline_pos

"""Get position on field line that passes r0 = (x0,y0,z0), at height z, for a given RD config c"""
function fieldline_pos(c::RotationalDiscontinuity, x0, y0, z0, z)
    int_sin(β, z) = sind(β * tanh(z))
    int_cos(β, z) = cosd(β * tanh(z))

    Ix = quadgk(ζ -> int_sin(c.β, ζ), z0, z)[1]
    Iy = quadgk(ζ -> int_cos(c.β, ζ), z0, z)[1]

    x = x0 + tand(c.θ) * Ix
    y = y0 + c.sign * tand(c.θ) * Iy
    return x, y, z
end

# This is not stable for large z. As tanh(z) -> 1, s₋(z) -> 0, and Ci(0) = -Inf
# Integrate[Sin[\[Beta] Tanh[z]], z]
# Integrate[Cos[\[Beta] Tanh[z]], z]
function _fieldline_pos(c::RotationalDiscontinuity, x0, y0, z0, z)
    s₊(z) = 1 + tanh(z)
    s₋(z) = 1 - tanh(z)
    function x_indefinite(z)
        β = deg2rad(c.β)
        part1 = sin(β) * (-cosint(β * s₋(z)) - cosint(β * s₊(z)))
        part2 = cos(β) * (sinint(β * s₋(z)) + sinint(β * s₊(z)))
        return 1 / 2 * (part1 + part2)
    end
    function y_indefinite(z)
        β = deg2rad(c.β)
        part1 = sin(β) * (sinint(β * s₊(z)) - sinint(β * s₋(z)))
        part2 = cos(β) * (cosint(β * s₊(z)) - cosint(β * s₋(z)))
        return 1 / 2 * (part1 + part2)
    end
    Ix = x_indefinite(z) - x_indefinite(z0)
    Iy = y_indefinite(z) - y_indefinite(z0)

    x = x0 + tand(c.θ) * Ix
    y = y0 + c.sign * tand(c.θ) * Iy
    return x, y, z
end

fieldline_pos(c, r0, z) = fieldline_pos(c, r0[1], r0[2], r0[3], z)

# SpecialFunctions
