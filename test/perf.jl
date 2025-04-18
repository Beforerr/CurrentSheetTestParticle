using CurrentSheetTestParticle
using BenchmarkTools

save_everystep = false
verbose = false
dtmax = 1e-2
diffeq = (; verbose, dtmax)

p = RDProblemParams(θ=45, β=90, v=8)

@b solve_params($p)

@benchmark solve_params($p; f=trace_normalized_B!, diffeq...) seconds = 1
@benchmark solve_params($p; f=trace_normalized_B, diffeq...) seconds = 1