using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin

    # wave parameters
    d, H, L, g, rho = 1, 0.1, 1, G, RHO # depth, height, length, gravity acceleration
    N = 40
    # Test: no mass transport in a flume
    k = 2π / L
    @time w1, df = fourier_approx(d, H, L; cc=2, N=N)

    @test w1.C ≈ w1.U # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / w1.C * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2, N=N)
    
    @test w1.C ≈ w.C

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
    @test 0 ≈ w.v.z(0,      w.eta.point(0))
    
    #velocity at wave lowest point
    @test 1e-15 > abs(w.v.z(π,w.eta.point(N)))

    #velocity inside a wave
    @test 1e-4 < abs(w.v.z(π/N,   w.eta.point(1)/2))

    #Test: vertical_velocity
    k = Output.wave_number(w, df)
    @test 0 ≈ vertical_velocity(w, 0,0, df)

    X = π/N
    Z = w.eta.point(1)/2
    @test vertical_velocity(w, X/k, Z/k, df) ≈ sqrt(g/k) * w.v.z(X, Z)
    
    #Test: horizontal_velocity
    @test horizontal_velocity(w, X/k, Z/k, df) ≈ sqrt(g/k) * w.v.x(X, Z)

    #Test: pressure
    @test pressure(w, X/k, Z/k, df) ≈ rho * g / k * pressure(w, X, Z)

    @test 1e-12 > abs(pressure(w,0,w.eta.point(0)))

    # Test: dispersion relation
    ω = √(g * k * tanh(k * d))
    @test k ≈ linear_wave_number(d, ω)

    # Test: elevation
    @test abs(elevation(w, 0.0)-w.eta.point(0)) < 1e-10

    for i in 1:N-1
        @test abs(elevation(w, pi*i/N)-w.eta.point(i)) < 1e-10
    end

    @test abs(elevation(w, pi)-w.eta.point(N)) < 1e-10

end
