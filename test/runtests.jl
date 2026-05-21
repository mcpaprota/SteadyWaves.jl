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


    @test SteadyWaves.Indirect.indirect_surface_tension(1,1,0) ≈ 0

    t = 1
    u = 1/t

    w,_ = SteadyWaves.linear_solution(1,0.1,1)

    @test Output.elevation(w,-1) ≈ Output.elevation(w,0,t,u)

    @test Output.elevation(w,-w.C) ≈ Output.elevation(w,0,t)


    @test Output.surface_tension(w,-1) ≈ Output.surface_tension(w,0,t,u)

    @test Output.surface_tension(w,-w.C) ≈ Output.surface_tension(w,0,t)


    @test Output.horizontal_velocity(w,-1,1) ≈ (Output.horizontal_velocity(w,0,1,t,u) -u)

    @test (Output.horizontal_velocity(w,0,1,t/w.C) -w.C) ≈ (Output.horizontal_velocity(w,0,1,t,u) -u)


    @test Output.vertical_velocity(w,-1,1) ≈ Output.vertical_velocity(w,0,1,t,u)

    @test Output.vertical_velocity(w,0,1,t/w.C) ≈ Output.vertical_velocity(w,0,1,t,u)

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

    @test Linear.test_solution_for_linearity(w1) > 1e-5

    @test w1.C ≈ w1.U # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / w1.C * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2, N=N,
        eta_type = SteadyWaves.Params.DIRECT_ELEVATION
    )

    @test Linear.test_solution_for_linearity(w1) > 1e-5
    
    @test w1.C ≈ w.C

    wd = dimensional(w,df)

    # Test: wave_length
    @test L ≈ wavelength(wd)

    # Test: wave_period
    @test T ≈ wave_period(wd)

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
    k = Output.wave_number(wd)
    @test 0 ≈ vertical_velocity(wd, 0,0)

    X = π/N
    Z = w.eta.point(1)/2
    @test vertical_velocity(wd, X/k, Z/k) ≈ sqrt(g/k) * w.v.z(X, Z)
    
    #Test: horizontal_velocity
    @test horizontal_velocity(wd, X/k, Z/k) ≈ sqrt(g/k) * w.v.x(X, Z)

    #Test: pressure
    @test pressure(wd, X/k, Z/k) ≈ rho * g / k * pressure(w, X, Z)

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

    @test SteadyWaves.Output.surface_tension(w,0) < 0

    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=SteadyWaves.Params.CC_EULER, N=N,
        eta_type = SteadyWaves.Params.DIRECT_ELEVATION,
        deep_water = true

    )
    w2,df2 = fourier_approx(1, 0.1, L; pc=2, cc=2,N = N, 
    eta_type=SteadyWaves.Params.DIRECT_ELEVATION,
    wave_type=SteadyWaves.Params.GRAVITY_CAPILLARY_WAVE
    )


    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=SteadyWaves.Params.CC_EULER, N=N,
        eta_type = SteadyWaves.Params.DIRECT_ELEVATION,
        deep_water = Params.STABLE
    )
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

    @test Linear.test_solution_for_linearity(w1) > 1e-5

    @test w1.C ≈ w1.U # c√(k/g) = Ū√(k/g)

    # Test: c = L/T
    T = L / w1.C * √(k / 9.81) # T = c/L
    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=2, N=N,
        eta_type = SteadyWaves.Params.FOURIER_ELEVATION
    )
    
    @test Linear.test_solution_for_linearity(w1) > 1e-5

    @test w1.C ≈ w.C
    wd = dimensional(w,df)

    # Test: wave_length
    @test L ≈ wavelength(wd)

    # Test: wave_period
    @test T ≈ wave_period(wd)

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
    k = Output.wave_number(wd)
    @test 0 ≈ vertical_velocity(wd, 0,0)

    X = π/N
    Z = w.eta.point(1)/2
    @test vertical_velocity(wd, X/k, Z/k) ≈ sqrt(g/k) * w.v.z(X, Z)
    
    #Test: horizontal_velocity
    @test horizontal_velocity(wd, X/k, Z/k) ≈ sqrt(g/k) * w.v.x(X, Z)

    #Test: pressure
    @test pressure(wd, X/k, Z/k) ≈ rho * g / k * pressure(w, X, Z)

    @test 1e-12 > abs(pressure(w,0,w.eta.point(0)))

    @test SteadyWaves.Output.surface_tension(w,0) < 0

    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=SteadyWaves.Params.CC_EULER, N=N,
        eta_type = SteadyWaves.Params.FOURIER_ELEVATION,
        deep_water = true

    )

    @time w, df = fourier_approx(1, 0.1, T; pc=2, cc=SteadyWaves.Params.CC_EULER, N=N,
        eta_type = SteadyWaves.Params.FOURIER_ELEVATION,
        deep_water = Params.STABLE
    )
end

@testset "SteadyWaves.jl - still_water" begin
    
    d = 1
    L = 1
    H = 0.1
    N = 40
    
    @time ws, df = Steady.fourier_approx(d,H,L,nothing,Params.ConfigStruct(
        eta_type=Params.FOURIER_ELEVATION,
        cc=Params.CC_STOKES,
        pc=Params.PC_LENGTH,
        reference_level=Params.STILL_DEPTH
        ),
        Physics.DEFAULT_PHYSICS,
        N=N
    )
    
    ws = SteadyWaves.dimensional(ws,df)
    print(ws.eta.avg)

    @time wm, df = Steady.fourier_approx(ws.eta.avg,H,L,nothing,Params.ConfigStruct(
        eta_type=Params.FOURIER_ELEVATION,
        cc=Params.CC_STOKES,
        pc=Params.PC_LENGTH,
        reference_level=Params.MEAN_DEPTH
        ),
        Physics.DEFAULT_PHYSICS,
        N=N
    )

    wm = SteadyWaves.dimensional(wm,df)

    @test sum(abs.(ws.eta.a .- wm.eta.a)) < 1e-10

    @test sum(abs.(ws.v.b .- wm.v.b)) < 1e-10

    @test ws.C ≈ wm.C
    
    @test ws.Q ≈ wm.Q

    @test ws.R ≈ wm.R

    @test ws.U ≈ wm.U
end


@testset "SteadyWaves.jl - fourier capillary" begin
    
    N = 20

    L = 0.01

    H = 0.2*L

    eta_type = SteadyWaves.Params.FOURIER_ELEVATION

    wg,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type= eta_type,
    wave_type=SteadyWaves.Params.GRAVITY_WAVE,
    deep_water = true
    )

    wc,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type=eta_type,
    wave_type=SteadyWaves.Params.CAPILLARY_WAVE,
    deep_water = true
    )

    wgc,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type=eta_type,
    wave_type=SteadyWaves.Params.GRAVITY_CAPILLARY_WAVE,
    deep_water = true
    )

    @test Linear.test_solution_for_linearity(wg) > 1e-5

    @test Linear.test_solution_for_linearity(wc) > 1e-5

    @test Linear.test_solution_for_linearity(wgc) > 1e-5

    @test abs(Crapper.crapper_test_of_solution(wgc)) < abs(Crapper.crapper_test_of_solution(wg))

    @test abs(Crapper.crapper_test_of_solution(wc)) < abs(Crapper.crapper_test_of_solution(wgc))

end


@testset "SteadyWaves.jl - direct capillary" begin
    
    N = 20

    L = 0.01

    H = 0.2*L

    eta_type = SteadyWaves.Params.DIRECT_ELEVATION

    wg,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type= eta_type,
    wave_type=SteadyWaves.Params.GRAVITY_WAVE,
    deep_water = true
    )

    wc,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type=eta_type,
    wave_type=SteadyWaves.Params.CAPILLARY_WAVE,
    deep_water = true
    )

    wgc,_ = fourier_approx(1, H, L; pc=1, cc=2,N = N, 
    eta_type=eta_type,
    wave_type=SteadyWaves.Params.GRAVITY_CAPILLARY_WAVE,
    deep_water = true
    )

    @test Linear.test_solution_for_linearity(wg) > 1e-5

    @test Linear.test_solution_for_linearity(wc) > 1e-5

    @test Linear.test_solution_for_linearity(wgc) > 1e-5

    @test abs(Crapper.crapper_test_of_solution(wgc)) < abs(Crapper.crapper_test_of_solution(wg))

    @test abs(Crapper.crapper_test_of_solution(wc)) < abs(Crapper.crapper_test_of_solution(wgc))

end

@testset "SteadyWaves.jl - dimensionless input" begin
    N = 32
    kd = 1
    kH = 0.2

    k = 1
    L = 2pi/k

    d = kd/k
    H = kH/k


    @time w, _ = Steady.fourier_approx(d,H,L,nothing,Params.ConfigStruct(
        eta_type=Params.FOURIER_ELEVATION,
        cc=Params.CC_EULER,
        pc=Params.PC_LENGTH,
        ),
        Physics.DEFAULT_PHYSICS,
        N=N
    )

    @time wd, _ = Steady.dimensionless_fourier_approx(kd,kH,Params.ConfigStruct(
        eta_type=Params.FOURIER_ELEVATION,
        cc=Params.CC_EULER
        ),
        dimensionless_sigma = 0,
        N = N
    )

    @test sum(abs.(w.eta.a .- wd.eta.a)) < 1e-10

    @test w.D ≈ wd.D

    @test w.C ≈ wd.C
    
    @test w.Q ≈ wd.Q

    @test w.R ≈ wd.R

    @test w.U ≈ wd.U

    @time wd, _ = Steady.dimensionless_fourier_approx(kd,kH,Params.ConfigStruct(
        eta_type=Params.FOURIER_ELEVATION,
        cc=Params.CC_EULER,
        indirect_celerity = true,
        ),
        dimensionless_sigma = 0,
        N = N
    )

    @test sum(abs.(w.eta.a .- wd.eta.a)) < 1e-10

    @test w.D ≈ wd.D

    @test w.C ≈ wd.C
    
    @test w.Q ≈ wd.Q

    @test w.R ≈ wd.R

    @test w.U ≈ wd.U

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