F0func(x) = SVector(0.0, 0.0, 0.0)

@kwdef struct Parameters{TE,TB,TF,Q,M}
    E::TE = F0func
    B::TB = F0func
    F::TF = F0func
    q::Q = 1.0
    m::M = 1.0
end

Bfunc(t::Tuple) = t[3]
Bfunc(p::Parameters) = p.B