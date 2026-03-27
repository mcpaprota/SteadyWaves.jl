module Wave

using ..Index:IndexStruct
using ..Surface:SurfaceStruct
using ..Velocity:VelocityStruct

struct WaveStruct
    eta
    v
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
    P
    raw # u array to keep compatibility with not updated functions (u,N instead of w)
    WaveStruct(eta,v,D,C,R,H,U,Q,N,L,T,F,P,raw) = new(eta,v,D,C,R,H,U,Q,N,L,T,F,P,raw)

    # create struct that replaces values given as key words in default
    WaveStruct(
        default::WaveStruct=nothing;
        eta = nothing,
        v = nothing,
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
        P = nothing,
        raw = nothing,
    ) = new(
        eta === nothing ? default.eta : eta,
        v === nothing ? default.v : v,
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
        P === nothing ? default.P : P,
        raw === nothing ? default.raw : raw,
    )

    WaveStruct(
        default::WaveStruct,
        df::WaveStruct;
        eta = nothing,
        v = nothing,
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
        P = nothing,
        raw = nothing,
    ) = new(
        eta === nothing ? default.eta : (w_c, u) -> eta * df.eta(w_c,u),
        v === nothing ? default.v : (w_c, u) -> v * df.v(w_c,u),
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
        P   === nothing ? default.P   : (w_c, u) -> P   * df.P(w_c,u),
        raw === nothing ? default.raw : (w_c, u) -> raw * df.raw(w_c,u),
    )

    # creates compiler X = (w_c, u) -> u[idx.X] from IndexStruct
    WaveStruct(idx::IndexStruct) = new(
        SurfaceStruct(idx),
        VelocityStruct(idx),
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
        (w_c, u) -> nothing,
        (w_c, u) -> u,
    )

    # create structure from array u and compiler X = compiler.X(compiler,u)
    WaveStruct(u,compiler::WaveStruct,inner_compiler::WaveStruct) = new(
        compiler.eta(inner_compiler,u),
        compiler.v(inner_compiler,u),
        compiler.D(inner_compiler,u),
        compiler.C(inner_compiler,u),
        compiler.R(inner_compiler,u),
        compiler.H(inner_compiler,u),
        compiler.U(inner_compiler,u),
        compiler.Q(inner_compiler,u),
        compiler.N(inner_compiler,u),
        compiler.L(inner_compiler,u),
        compiler.T(inner_compiler,u),
        compiler.F(inner_compiler,u),
        compiler.P(inner_compiler,u),
        compiler.raw(inner_compiler,u),
    )

    WaveStruct(u,compiler::WaveStruct) = WaveStruct(u,compiler,compiler)
end

end