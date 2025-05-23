module TestParticleAnalysis

# Write your package code here.
using AlgebraOfGraphics
using LaTeXStrings

tp_pairing() = (
    β=:β => L"β",
    μ0=:μ0 => L"μ_0",
    μ0_n=:μ0 => x -> string(round(x, digits=2)),
    t_cs=:t1 => L"T_{cross}",
    T=:T,
    Δs_perp=:Δs_perp => L"Δs_{\perp}",
    Δs_para=:Δs_para => L"Δs_{\parallel}",
    κ_para=:κ_para => L"κ_{\parallel}",
    κ_perp=:κ_perp => L"κ_{\perp}",
    κ_ratio=:κ_ratio => L"κ_{\perp}/κ_{\parallel}",
)

tp_mapping(m) = (
    t_cs_μ0=mapping(m.t_cs) * AlgebraOfGraphics.density(),
    Δt_μ0=mapping(:Δt) * AlgebraOfGraphics.density(),
    μ=mapping(:μ0, :μ1),
    s_perp=mapping(m.Δs_perp) * AlgebraOfGraphics.density(),
    s_para=mapping(m.Δs_para) * AlgebraOfGraphics.density(),
    κ_para=mapping(m.κ_para) * AlgebraOfGraphics.density(),
    κ_perp=mapping(m.κ_perp) * AlgebraOfGraphics.density(),
    κs=mapping([m.κ_perp, m.κ_para, m.κ_ratio], row=dims(1) => renamer(["", "", ""])),
    v_c=mapping(; color=:v => nonnumeric)
)




end
