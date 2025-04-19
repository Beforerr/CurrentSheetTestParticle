using ChairmarksForAirspeedVelocity
using CurrentSheetTestParticle
using Logging

const SUITE = BenchmarkGroup()

p = RDProblemParams(θ=45, β=90, v=8)
SUITE["main"]["RotationalDiscontinuity"] = @benchmarkable solve_params($p)
SUITE["main"]["RotationalDiscontinuity - No save"] = @benchmarkable solve_params($p, save_everystep=false)
SUITE["main"]["RotationalDiscontinuity - No save - NullLogger"] = @benchmarkable with_logger(NullLogger()) do
    solve_params($p, save_everystep=false)
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
