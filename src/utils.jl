# Accessors
𝓑(t::Tuple) = t[3]
𝓑(p::Parameters) = p.B
𝓑(f::Function) = f

_alg(x) = x
function _alg(alg::Symbol)
    return @match alg begin
        :AutoVern9 => AutoVern9(Rodas4P())
        :Boris => error("Boris algorithm is not implemented yet")
        _ => error("Unknown algorithm: $alg")
    end
end

function produce_or_load(f, path; load = true)
    return if load && isfile(path)
        JLD2.load(path, "sol")
    else
        sol = f()
        load && begin
            mkpath(dirname(path))
            JLD2.jldsave(path; sol)
        end
        sol
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

distance(A::AbstractMatrix, B::AbstractMatrix; dims = 1) = distance(eachslice(A, dims = dims), eachslice(B, dims = dims))

"""
Calculate the distance between two parallel lines.

Each line is defined by one point and a direction vector.
"""
function distance(p1, p2, d)
    r = p2 - p1
    # Calculate the value under the square root
    val = r ⋅ r - (r ⋅ d)^2 / (d ⋅ d)
    # Handle potential floating point errors that can cause slightly negative values
    return sqrt(max(0.0, val))
end

get_r(u) = u[SA[1, 2, 3]]
get_z(u) = u[3]
get_v(u) = u[SA[4, 5, 6]]
_get_r(u) = u[SA[1, 2, 3]]
get_q2m(p) = p[1] / p[2]
get_B(p) = p[3]

function Larmor_vector(𝐯, 𝐁, q2m)
    return (𝐁 × 𝐯) ./ (q2m * 𝐁 ⋅ 𝐁)
end

guiding_center(𝐫, 𝐯, 𝐁, q2m) = 𝐫 - Larmor_vector(𝐯, 𝐁, q2m)
guiding_center(xu, Bf::Function, q2m = 1) =
    guiding_center(get_r(xu), get_v(xu), Bf(xu), q2m)

guiding_centers(sol::ODESolution) = BroadcastArray(guiding_center, sol.u, 𝓑(sol))


energy(v; m = mp) = 0.5 * m * v^2 |> u"MeV"
velocity(E; m = mp) = sqrt(2 * E / m) |> upreferred
velocity(v::Velocity) = v
