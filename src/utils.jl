_alg(x) = x
function _alg(alg::Symbol)
    @match alg begin
        :AutoVern9 => AutoVern9(Rodas4P())
        :Boris => error("Boris algorithm is not implemented yet")
        _ => error("Unknown algorithm: $alg")
    end
end

function distance(r1s, r2s)
    result = minimum(r1s) do r1
        minimum(r2s) do r2
            (r1 - r2) ⋅ (r1 - r2)
        end
    end
    return sqrt(result)
end

distance(A::AbstractMatrix, B::AbstractMatrix) = distance(eachcol(A), eachcol(B))

"""
Calculate the distance between two parallel lines.

Each line is defined by one point and a direction vector.
"""
function distance(p1, p2, d)
    v = p2 .- p1
    n = cross(v, d)
    return norm(n) / norm(d)
end

get_r(u) = u[SA[1, 2, 3]]
get_z(u) = u[3]
get_v(u) = @view u[4:6]
_get_r(u) = u[SA[1, 2, 3]]
get_q2m(p) = p[1] / p[2]
get_B(p) = p[3]

function Larmor_vector(𝐯, 𝐁, q2m)
    (𝐁 × 𝐯) ./ (q2m * 𝐁 ⋅ 𝐁)
end

guiding_center(𝐫, 𝐯, 𝐁, q2m) = 𝐫 - Larmor_vector(𝐯, 𝐁, q2m)
guiding_center(xu, Bf::Function, q2m=1) =
    guiding_center(get_r(xu), get_v(xu), Bf(xu), q2m)