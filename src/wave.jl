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
        default::WaveStruct;
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

end