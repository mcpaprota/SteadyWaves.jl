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
        eta === nothing ? default.eta : (w_c, u) -> eta * df.eta(w_c,u),
        psi === nothing ? default.psi : (w_c, u) -> psi * df.psi(w_c,u),
        D   === nothing ? default.D   : (w_c, u) -> D   * df.D(w_c,u),
        C   === nothing ? default.C   : (w_c, u) -> C   * df.C(w_c,u),
        R   === nothing ? default.R   : (w_c, u) -> R   * df.R(w_c,u),
        H   === nothing ? default.H   : (w_c, u) -> H   * df.H(w_c,u),
        U   === nothing ? default.U   : (w_c, u) -> U   * df.U(w_c,u),
        Q   === nothing ? default.Q   : (w_c, u) -> Q   * df.Q(w_c,u),
        N   === nothing ? default.N   : (w_c, u) -> N   * df.N(w_c,u),
        L   === nothing ? default.L   : (w_c, u) -> L   * df.L(w_c,u),
        T   === nothing ? default.T   : (w_c, u) -> T   * df.T(w_c,u),
        F   === nothing ? default.F   : (w_c, u) -> F   * df.F(w_c,u),
        raw === nothing ? default.raw : (w_c, u) -> raw * df.raw(w_c,u),
    )

    # creates compiler X = (w_c, u) -> u[idx.X] from IndexStruct
    WaveStruct(idx::IndexStruct) = new(
        SurfaceStruct(idx),
        (w_c, u) -> u[idx.psi],
        (w_c, u) -> u[idx.D],
        (w_c, u) -> u[idx.C],
        (w_c, u) -> u[idx.R],
        (w_c, u) -> u[idx.H],
        (w_c, u) -> u[idx.U],
        (w_c, u) -> u[idx.Q],
        (w_c, u) -> idx.N,
        (w_c, u) -> nothing,
        (w_c, u) -> nothing,
        (w_c, u) -> nothing,
        (w_c, u) -> u,
    )

    # create structure from array u and compiler X = compiler.X(compiler,u)
    WaveStruct(u,compiler::WaveStruct) = new(
        compiler.eta(compiler,u),
        compiler.psi(compiler,u),
        compiler.D(compiler,u),
        compiler.C(compiler,u),
        compiler.R(compiler,u),
        compiler.H(compiler,u),
        compiler.U(compiler,u),
        compiler.Q(compiler,u),
        compiler.N(compiler,u),
        compiler.L(compiler,u),
        compiler.T(compiler,u),
        compiler.F(compiler,u),
        compiler.raw(compiler,u),
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
function dimensional_factor_compiler(d,g,rho) 

    return WaveStruct(
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # eta
	    (w_c, u) -> bernoulli_factor(  w_c.D(w_c,u), d, g),	    # psi
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # D
	    (w_c, u) -> speed_factor(      w_c.D(w_c,u), d, g),	    # C
	    (w_c, u) -> bernoulli_factor(  w_c.D(w_c,u), d, g),	    # R
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # H
	    (w_c, u) -> speed_factor(      w_c.D(w_c,u), d, g),	    # U
	    (w_c, u) -> flux_factor(       w_c.D(w_c,u), d, g),	    # Q
	    (w_c, u) -> 1,	                                        # N
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # L
	    (w_c, u) -> period_factor(     w_c.D(w_c,u), d, g),	    # T
	    (w_c, u) -> power_factor(      w_c.D(w_c,u), d, g, rho),	# F
	    (w_c, u) -> 1	                                        # raw
    )
end

end