using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin

    # wave parameters
    d, H, L, g, rho = 1, 0.1, 1, G, RHO # depth, height, length, gravity acceleration
    N = 10
    # Test: no mass transport in a flume
    k = 2π / L
    @time w, df = fourier_approx(d, H, L; cc=2)
    u1 = w.raw

    @test u1[22] ≈ u1[26] # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / u1[22] * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2)
    u2 = w.raw
    
    @test u1[22] ≈ u2[22]

    # Test: wave_length
    @test L ≈ wavelength(w, df)

    # Test: wave_period
    @test T ≈ wave_period(w, df)

    # Test: shoaling
    K = topo_approx([d, d], H, L)
    @test K[1] ≈ K[2]
 # Test: w.v.z
    #velocity at 3 bottom points
    @test 0 ≈ w.v.z(0,     0)
    @test 0 ≈ w.v.z(π/2,   0)
    @test 0 ≈ w.v.z(π,     0)

    #velocity at wave highest point
    @test 0 ≈ w.v.z(0,      u2[1])
    
    #velocity at wave lowest point
    @test 1e-15 > abs(w.v.z(π,     u2[N+1]))

    #velocity inside a wave
    @test 1e-4 < abs(w.v.z(π/N,   u2[2]/2))

    #Test: vertical_velocity
    k = Output.wave_number(w, df)
    @test 0 ≈ vertical_velocity(w, 0,0, df)

    X = π/N
    Z = u2[2]/2
    @test vertical_velocity(w, X/k, Z/k, df) ≈ sqrt(g/k) * w.v.z(X, Z)
    
    #Test: horizontal_velocity
    @test horizontal_velocity(w, X/k, Z/k, df) ≈ sqrt(g/k) * w.v.x(X, Z)

    #Test: pressure
    @test pressure(w, X/k, Z/k, df) ≈ rho * g / k * pressure(w, X, Z)

    @test 1e-12 > abs(pressure(w,0,u2[1]))

    # Test: dispersion relation
    ω = √(g * k * tanh(k * d))
    @test k ≈ linear_wave_number(d, ω)
    
    w2,df2 = fourier_approx(1, 0.1, T; pc=2, cc=2, sigma=SteadyWaves.Physics.SIGMA)

    @test w2.eta.point(0) /df2.H < w.eta.point(0) /df.H

    @test w2.eta.point(N) /df2.H < w.eta.point(N) /df.H

    @test w2.eta.point(Int(N/2)) /df2.H > w.eta.point(Int(N/2)) /df.H
end
