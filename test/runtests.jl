using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin
    # no mass transport in a flume
    d, H, L, g = 1, 0.1, 1, 9.81 # depth, height, length, gravity acceleration
    k = 2π / L
    u1 = fourier_approx(d, H, L; cc=2)
    @test u1[22] ≈ u1[26] # c√(k/g) = Ū√(k/g)
    # c = L/T
    T = L / u1[22] * √(k / g) # T = c/L
    u2 = fourier_approx(1, 0.1, T; pc=2, cc=2)
    @test u1[22] ≈ u2[22]
    # dispersion relation
    ω = √(g * k * tanh(k * d))
    @test k ≈ SteadyWaves.dispertion_relation(d, ω)
end
