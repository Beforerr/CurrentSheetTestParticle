const DEFAULT_θ = 45
const DEFAULT_β = 90
const DEFAULT_SIGN = 1

abstract type AbstractNormalization end

struct Normalization{L, V, T} <: AbstractNormalization
    L₀::L
    V₀::V
    T₀::T
end


"""
# Keyword Arguments
- `dir` : direction where the field depends on, 1, 2, or 3
"""
B(r, conf::RotationalDiscontinuity; dir = 3) = conf(r; dir)

function normalization(Π::RotationalDiscontinuity; m = mp)
    L₀ = Π.L
    Ω₀ = q * Π.B / m
    T₀ = 1 / Ω₀ |> upreferred
    V₀ = L₀ / T₀ |> upreferred
    return Normalization(L₀, V₀, T₀)
end
