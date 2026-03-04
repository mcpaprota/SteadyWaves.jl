using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin

    # wave parameters
    d, H, L, g = 1, 0.1, 1, G # depth, height, length, gravity acceleration

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
    @test k ≈linear_wave_number(d, ω)

    # Test: shoaling
    K = topo_approx([d, d], H, L)
    @test K[1] ≈ K[2]

    # Test: dimensionless_vertical_velocity
    #velocity at 3 bottom points
    @test 0 ≈ dimensionless_vertical_velocity(u2,N,     0,     0)
    @test 0 ≈ dimensionless_vertical_velocity(u2,N,     π/2,   0)
    @test 0 ≈ dimensionless_vertical_velocity(u2,N,     π,     0)

    #velocity at wave highest point
    @test 0 ≈ dimensionless_vertical_velocity(u2, N,    0,      u2[1])
    
    #velocity at wave lowest point
    @test 1e-15 > abs(dimensionless_vertical_velocity(u2, N,     π,     u2[N+1]))

    #velocity inside a wave
    @test 1e-4 < abs(dimensionless_vertical_velocity(u2, N,   π/N,   u2[2]/2))

    #Test: vertical_velocity
    k = Output.wave_number(u2,d,N)
    @test 0 ≈ vertical_velocity(u2,N,0,0,k)

    X = π/N
    Z = u2[2]/2
    @test vertical_velocity(u2,N, X/k, Z/k,k) ≈ sqrt(k/g) * dimensionless_vertical_velocity(u2, N, X, Z)

    #Test: horizontal_velocity
    @test horizontal_velocity(u2,N, X/k, Z/k,k) ≈ sqrt(k/g) * dimensionless_horizontal_velocity(u2, N, X, Z)

 
end
