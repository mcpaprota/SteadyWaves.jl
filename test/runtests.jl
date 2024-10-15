using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin

    # wave parameters
    d, H, L, g = 1, 0.1, 1, 9.81 # depth, height, length, gravity acceleration

    # Test: no mass transport in a flume
    k = 2π / L
    u1 = fourier_approx(d, H, L; cc=2)
    @test u1[22] ≈ u1[26] # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / u1[22] * √(k / g) # T = c/L
    u2 = fourier_approx(1, 0.1, T; pc=2, cc=2)
    @test u1[22] ≈ u2[22]

    # Test: wave_length
    @test L ≈ wavelength(u2, d, 10)

    # Test: wave_period
    @test T ≈ wave_period(u2, d, 10)

    # Test: dispersion relation
    ω = √(g * k * tanh(k * d))
    @test k ≈ dispersion_relation(d, ω)

    # Test: shoaling
    K = shoaling_approx([d, d], H, L)
    @test K[1] ≈ K[2]

end
