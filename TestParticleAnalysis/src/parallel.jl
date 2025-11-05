module ParallelTransport

using LinearAlgebra
using DimensionalData
using FiniteDifferences
using InverseLaplace

export stationary_distribution, msd

unwrap(x) = parent(x)
unwrap(x::Array) = x

"""
    stationary_distribution(P)

Compute stationary distribution π for column-stochastic P (P * π = π).
Returns a real, normalized vector.
"""
function stationary_distribution(P::AbstractMatrix)
    sum_close_to_one(x) = isapprox(sum(x), 1)
    @assert all(sum_close_to_one, eachcol(P))
    vals, vecs = eigen(unwrap(P))
    k = argmin(abs.(vals .- 1))
    π = real.(vecs[:, k])
    π ./= sum(π)
    return π
end

# uniform in magnetic moment ⇒ pitch-angle cosine PDF ∝ |μ|,
stationary_distribution(μs::AbstractVector) = abs.(μs)

"""
    A_matrix(P, v, Psi, k, s)

Build A(k,s) with elements
A[i,j] = P[i,j] * Ψ_j(s - i k v_j),
where Ψ_j is the Laplace transform of ψ_j(τ).
- P: M×M, column-stochastic (from j → to i).
- v: length-M vector of velocities v_j.
- Psi: Vector of functions Psi[j](λ) = Ψ_j(λ).
"""
function A_matrix(P, v, Psi, k, s)
    M = size(P, 1)
    A = zeros(ComplexF64, M, M)
    μs = P.dims[1]
    for j in 1:M  # start state
        μ_j = μs[j]
        λ = s - im * k * v * μ_j
        Ψ = Psi(λ, μ_j)
        for i in 1:M  # end state
            A[i, j] = P[i, j] * Ψ
        end
    end
    return A
end

Φ̃(λ, μ, Ψ̃) = (1 - Ψ̃(λ, μ)) / λ

"""
    D_matrix(v, Psi, k, s)

Build diagonal D(k,s) with entries
D[i,i] = Φ̃_i(s - i k v_i) = (1 - Ψ_i(s - i k v_i)) / (s - i k v_i).
"""
function D_matrix(μs, v, Ψ̃, k, s)
    M = length(μs)
    D = zeros(ComplexF64, M, M)
    for i in 1:M
        μ = μs[i]
        λ = s - im * k * v * μ
        D[i, i] = Φ̃(λ, μ, Ψ̃)
    end
    return D
end

"""
    F_tilde(k, s, P, v, Ψ̃; π = nothing)

Compute F̂(k,s): Fourier-Laplace transform of total position PDF.
"""
function F_tilde(k, s, P, v, Ψ̃; π = nothing)
    π = @something π stationary_distribution(P)
    A = A_matrix(P, v, Ψ̃, k, s)
    μs = P.dims[1]
    D = D_matrix(μs, v, Ψ̃, k, s)
    J = (I - A) \ π
    return sum(D * J)
end

"""
    M2_tilde(s, P, v, Psi; π = nothing)

Compute Laplace transform of MSD, M̃₂(s),
using a symmetric finite-difference approximation:
M̃₂(s) = -∂²F̂/∂k² |_{k=0}.

`p`: Number of grid points used for the finite-difference approximation.
"""
function M2_tilde(s, P, v, Psi; π = nothing, p = 5)
    return -central_fdm(p, 2)(k -> F_tilde(k, s, P, v, Psi; π = π), 0)
end

# ================= Inverse Laplace (Gaver-Stehfest) ====================
"""
    msd(t, P, v, Psi; N=18)

Compute ⟨Δz²(t)⟩ for time `t`.
using M₂̃(s) and Stehfest inverse Laplace.
"""
function msd(t, P, v, Psi; π = nothing, N = 18)
    μs = P.dims[1].val
    π = @something π stationary_distribution(μs)
    F̃M2 = s -> M2_tilde(s, P, v, Psi; π = π)
    tf = Float64(t)
    # return ILT(F̃M2)(t)
    return gaverstehfest(F̃M2, tf)
    # return Weeks(F̃M2)(tf)
end

end # module
