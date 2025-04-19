using ChairmarksForAirspeedVelocity
using CurrentSheetTestParticle
using Logging
using ChairmarksForAirspeedVelocity: Chairmarks
Chairmarks.DEFAULTS.seconds = 1

const SUITE = BenchmarkGroup()

p = RDProblemParams(θ=45, β=90, v=8)
SUITE["main"]["RotationalDiscontinuity"] = @benchmarkable solve_params($p)
SUITE["main"]["RotationalDiscontinuity - No save"] = @benchmarkable solve_params($p, save_everystep=false)

SUITE["main"]["RotationalDiscontinuity - No save - NullLogger"] = with_logger(NullLogger()) do
    @benchmarkable solve_params($p, save_everystep=false)
end

@static if Sys.isapple()
    # Not working
    # Note: need DiffEqGPU#master
    # ERROR: Metal does not support Float64 values, try using Float32 instead
    using Metal
    backend = Metal.MetalBackend()
    with_logger(NullLogger()) do
        @be solve_params($p, save_everystep=false, ensemblealg=EnsembleGPUArray(backend))
    end
end
