using SteadyWaves
using Test

@testset "SteadyWaves.jl - helpers" begin
    # wave parameters
    d, H, L, g, rho = 1, 0.1, 1, G, RHO # depth, height, length, gravity acceleration
    # Test: no mass transport in a flume
    k = 2π / L

    # Test: dispersion relation
    ω = √(g * k * tanh(k * d))
    @test k ≈ linear_wave_number(d, ω)
end

@testset "SteadyWaves.jl - direct elevation" begin

    # wave parameters
    d, H, L, g, rho = 1, 0.1, 1, G, RHO # depth, height, length, gravity acceleration
    N = 40
    # Test: no mass transport in a flume
    k = 2π / L
    @time w1, df1 = fourier_approx(d, H, L; cc=2, N=N,
        eta_type = SteadyWaves.Params.DIRECT_ELEVATION
    )

    @test w1.C ≈ w1.U # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / w1.C * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2, N=N,
        eta_type = SteadyWaves.Params.DIRECT_ELEVATION
    )
    
    @test w1.C ≈ w.C

    # Test: wave_length
    @test L ≈ wavelength(w, df)

    # Test: wave_period
    @test T ≈ wave_period(w, df)

    # Test: shoaling
    K = topo_approx([d, d], H, L, eta_type=SteadyWaves.Params.DIRECT_ELEVATION)
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

    # Test: elevation
    @test abs(elevation(w, 0.0)-w.eta.point(0)) < 1e-10

    for i in 1:N-1
        @test abs(elevation(w, pi*i/N)-w.eta.point(i)) < 1e-10
    end

    @test abs(elevation(w, pi)-w.eta.point(N)) < 1e-10

    #Test: derivative
    @test w.eta.point.dz_dx_1(0) ≈ 0

    @test abs(w.eta.point.dz_dx_1(N)) < 1e-15

    @test w.eta.z.dz_dx_1(0) ≈ 0
    
    @test abs(w.eta.z.dz_dx_1(π)) < 1e-15

end

@testset "SteadyWaves.jl - fourier elevation" begin

    # wave parameters
    d, H, L, g, rho = 1, 0.1, 1, G, RHO # depth, height, length, gravity acceleration
    N = 40
    # Test: no mass transport in a flume
    k = 2π / L
    @time w1, df1 = fourier_approx(d, H, L; cc=2, N=N,
        eta_type = SteadyWaves.Params.FOURIER_ELEVATION
    )

    @test w1.C ≈ w1.U # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / w1.C * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2, N=N,
        eta_type = SteadyWaves.Params.FOURIER_ELEVATION
    )
    
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

end

@testset "SteadyWaves.jl - dimensional factor" begin

    using SteadyWaves.DimensionalFactor: DimensionalFactor as DF

    d, L, g, rho = 1, 1, G, RHO

    kd = 2pi/L * d

    @test DF.dimensional_factor(kd,d,g,rho;L=1) ≈ DF.distance_factor(kd,d)

    @test DF.dimensional_factor(kd,d,g,rho;T=1) ≈ DF.period_factor(kd,d,g)

    @test DF.dimensional_factor(kd,d,g,rho;L=1,T=-1) ≈ DF.speed_factor(kd,d,g)

    @test DF.dimensional_factor(kd,d,g,rho;L=1,T=-3,M=1) ≈ DF.power_factor(kd,d,g,rho)

    @test DF.dimensional_factor(kd,d,g,rho;L=-1,T=-2,M=1) ≈ DF.pressure_factor(kd,d,g,rho)
end