using TestParticleAnalysis
using Test

@testset "TestParticleAnalysis.jl" begin
    # Write your tests here.
end


@testset "fieldline_pos" begin
    # 1) trivial case: θ = 0 => field is purely along z, so x,y must not change
    c = RotationalDiscontinuity(1.0, 0.0, 2.0, 1)  # B=1, θ=0°
    x0, y0, z0 = 3.0, -1.5, 0.0
    @test fieldline_pos(c, x0, y0, z0, 5.0) == (x0, y0, 5.0)

    # 2) symmetry check: going from z0->z1 or z1->z0 should give opposite x,y shifts
    c3 = RotationalDiscontinuity(1.0, 40.0, 0.7, 1)
    xA, yA, zA = 0.0, 0.0, -1.0
    xB, yB, zB = fieldline_pos(c3, xA, yA, zA, +1.0)
    # now start at the new point and go back to -1.0, we should get original (within quad tol)
    xC, yC, zC = fieldline_pos(c3, xB, yB, zB, -1.0)
    @test xC ≈ xA atol = 1.0e-12
    @test yC ≈ yA atol = 1.0e-12
    @test zC ≈ zA

    @testset "fieldline_pos asymptotic straight-line" begin
        # choose parameters where the asymptotic direction is nontrivial
        c = RotationalDiscontinuity(1.0, 35.0, 50, 1)
        # start high up, where tanh(z) ≈ 1
        x0, y0, z0 = 0.3, -0.4, 100.0
        z = 120.0
        Δz = z - z0
        # numerical value from the function
        x_num, y_num, z_num = fieldline_pos(c, x0, y0, z0, z)

        # asymptotic straight-line prediction
        x_exp = x0 + tand(c.θ) * sind(c.β) * Δz
        y_exp = y0 + c.sign * tand(c.θ) * cosd(c.β) * Δz
        @test x_num ≈ x_exp atol = 1.0e-12
        @test y_num ≈ y_exp atol = 1.0e-12
        @test z_num == z
    end
end
