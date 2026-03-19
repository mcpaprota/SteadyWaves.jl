module Index
"""
`u[eta_indexes(N)]`: free surface elevations *kη*
"""
function eta_indexes(N)
return 1:N+1
end

"""
`u[psi_indexes(N)]` : stream function coefficients *B*
"""
function psi_indexes(N)
return N+2:2N+1
end

"""
`u[2N+C_INDEX]`: wave celerity *c√(k/g)*
"""
const C_INDEX::Int = 2

"""
`u[2N+D_INDEX]`: mean water depth *kη̄*
"""
const D_INDEX::Int = 3

"""
`u[2N+Q_INDEX]`: volume flux due to waves *q√(k³/g)*
"""
const Q_INDEX::Int = 4

"""
`u[2N+R_INDEX]`: Bernoulli constant *rk/g*
"""
const R_INDEX::Int = 5

"""
`u[2N+U_INDEX]`: mean flow velocity *Ū√(k/g)*
"""
const U_INDEX::Int = 6

"""
`u[2N+H_INDEX]`: wave height *kH*
"""
const H_INDEX::Int = 7

struct IndexStruct
    eta
    psi
    D
    C
    R
    H
    U
    Q
    N
end

function default_indexes(N)
    return IndexStruct(
        eta_indexes(N),
        psi_indexes(N),
        2N+D_INDEX,
        2N+C_INDEX,
        2N+R_INDEX,
        2N+H_INDEX,
        2N+U_INDEX,
        2N+Q_INDEX,
        N
    )
end

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
    raw # u array to keep compatibility with not updated functions (u,N instead of w)
    WaveStruct(eta,psi,D,C,R,H,U,Q,N,raw) = new(eta,psi,D,C,R,H,U,Q,N,raw)

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
        raw === nothing ? default.raw : raw,
    )

    # creates compiler X = u -> u[idx.X] from IndexStruct
    WaveStruct(idx::IndexStruct) = new(
        u -> u[idx.eta],
        u -> u[idx.psi],
        u -> u[idx.D],
        u -> u[idx.C],
        u -> u[idx.R],
        u -> u[idx.H],
        u -> u[idx.U],
        u -> u[idx.Q],
        u -> idx.N,
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
        compiler.raw(u),
    )

end

export C_INDEX, D_INDEX, H_INDEX, Q_INDEX, R_INDEX, U_INDEX

export eta_indexes,psi_indexes

end