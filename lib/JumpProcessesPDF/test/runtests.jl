using JumpProcessesPDF
using Test
using Aqua

@testset "JumpProcessesPDF.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(JumpProcessesPDF)
    end
    # Write your tests here.
end
