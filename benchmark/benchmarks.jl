using ChairmarksForAirspeedVelocity
using CurrentSheetTestParticle

const SUITE = BenchmarkGroup()

p = RDProblemParams(θ=45, β=90, v=8)
SUITE["main"]["RotationalDiscontinuity"] = @benchmarkable solve_params($p)
SUITE["main"]["RotationalDiscontinuity - No save"] = @benchmarkable solve_params($p, save_everystep=false)
