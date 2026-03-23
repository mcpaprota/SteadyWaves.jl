module Wave

using ..Index:IndexStruct
using ..Surface:SurfaceStruct

struct WaveStruct
    eta
    psi
    D
    C
    R
    H
    U
    Q
    N
    L
    T
    F
    raw # u array to keep compatibility with not updated functions (u,N instead of w)
    WaveStruct(eta,psi,D,C,R,H,U,Q,N,L,T,F,raw) = new(eta,psi,D,C,R,H,U,Q,N,L,T,F,raw)

    # create struct that replaces values given as key words in default
    WaveStruct(
        default::WaveStruct=nothing;
        eta = nothing,
        psi = nothing,
        D = nothing,
        C = nothing,
        R = nothing,
        H = nothing,
        U = nothing,
        Q = nothing,
        N = nothing,
        L = nothing,
        T = nothing,
        F = nothing,
        raw = nothing,
    ) = new(
        eta === nothing ? default.eta : eta,
        psi === nothing ? default.psi : psi,
        D === nothing ? default.D : D,
        C === nothing ? default.C : C,
        R === nothing ? default.R : R,
        H === nothing ? default.H : H,
        U === nothing ? default.U : U,
        Q === nothing ? default.Q : Q,
        N === nothing ? default.N : N,
        L === nothing ? default.L : L,
        T === nothing ? default.T : T,
        F === nothing ? default.F : F,
        raw === nothing ? default.raw : raw,
    )

        WaveStruct(
        default::WaveStruct,
        df::WaveStruct;
        eta = nothing,
        psi = nothing,
        D = nothing,
        C = nothing,
        R = nothing,
        H = nothing,
        U = nothing,
        Q = nothing,
        N = nothing,
        L = nothing,
        T = nothing,
        F = nothing,
        raw = nothing,
    ) = new(
        eta === nothing ? default.eta : u -> eta * df.eta(u),
        psi === nothing ? default.psi : u -> psi * df.psi(u),
        D   === nothing ? default.D   : u -> D   * df.D(u),
        C   === nothing ? default.C   : u -> C   * df.C(u),
        R   === nothing ? default.R   : u -> R   * df.R(u),
        H   === nothing ? default.H   : u -> H   * df.H(u),
        U   === nothing ? default.U   : u -> U   * df.U(u),
        Q   === nothing ? default.Q   : u -> Q   * df.Q(u),
        N   === nothing ? default.N   : u -> N   * df.N(u),
        L   === nothing ? default.L   : u -> L   * df.L(u),
        T   === nothing ? default.T   : u -> T   * df.T(u),
        F   === nothing ? default.F   : u -> F   * df.F(u),
        raw === nothing ? default.raw : u -> raw * df.raw(u),
    )

    # creates compiler X = u -> u[idx.X] from IndexStruct
    WaveStruct(idx::IndexStruct) = new(
        SurfaceStruct(idx),
        u -> u[idx.psi],
        u -> u[idx.D],
        u -> u[idx.C],
        u -> u[idx.R],
        u -> u[idx.H],
        u -> u[idx.U],
        u -> u[idx.Q],
        u -> idx.N,
        u -> nothing,
        u -> nothing,
        u -> nothing,
        u -> u,
    )

    # create structure from array u and compiler X = compiler.X(u)
    WaveStruct(u,compiler::WaveStruct) = new(
        compiler.eta(u),
        compiler.psi(u),
        compiler.D(u),
        compiler.C(u),
        compiler.R(u),
        compiler.H(u),
        compiler.U(u),
        compiler.Q(u),
        compiler.N(u),
        compiler.L(u),
        compiler.T(u),
        compiler.F(u),
        compiler.raw(u),
    )

end

function distance_factor(kd,d) 
    return kd/d
end

function period_factor(kd,d,g) 
    return sqrt(g * kd/d)
end

function speed_factor(kd,d,g)
    return sqrt(g * kd/d)
end

function power_factor(kd,d,g,rho)
    return rho * kd/d * sqrt((g * kd/d)^3)    
end

function bernoulli_factor(kd,d,g)
    return kd/d / g
end

function flux_factor(kd,d,g)
    return kd/d * sqrt(kd/d/ g) 
end

# returns compiler that produce factor to multiply dimensional values into dimentionless 
function dimensional_factor_compiler(u_to_kd,d,g,rho) 
    return WaveStruct(
	    u -> distance_factor(   u_to_kd(u), d),	        # eta
	    u -> bernoulli_factor(  u_to_kd(u), d, g),	    # psi
	    u -> distance_factor(   u_to_kd(u), d),	        # D
	    u -> speed_factor(      u_to_kd(u), d, g),	    # C
	    u -> bernoulli_factor(  u_to_kd(u), d, g),	    # R
	    u -> distance_factor(   u_to_kd(u), d),	        # H
	    u -> speed_factor(      u_to_kd(u), d, g),	    # U
	    u -> flux_factor(       u_to_kd(u), d, g),	    # Q
	    u -> 1,	                                        # N
	    u -> distance_factor(   u_to_kd(u), d),	        # L
	    u -> period_factor(     u_to_kd(u), d, g),	    # T
	    u -> power_factor(      u_to_kd(u), d, g, rho),	# F
	    u -> 1	                                        # raw
    )
end


end