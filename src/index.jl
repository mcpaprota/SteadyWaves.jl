module Index
"""
`u[elevation_indexes(N)]`: free surface elevations *kη*
"""
function elevation_indexes(N)
return 1:N+1
end

"""
`u[stream_indexes(N)]` : stream function coefficients *B*
"""
function stream_indexes(N)
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
    D
    C
    R
    H
    U
    Q
end

const INDEX_STRUCT = IndexStruct(
    N -> 2N+D_INDEX,
    N -> 2N+C_INDEX,
    N -> 2N+R_INDEX,
    N -> 2N+H_INDEX,
    N -> 2N+U_INDEX,
    N -> 2N+Q_INDEX,
)

struct InterpreterStruct
    eta
    D
    C
    R
    H
    U
    Q
    T
    F
    L

    InterpreterStruct(
        eta,
        id=INDEX_STRUCT;
        D=nothing,
        C=nothing,
        R=nothing,
        H=nothing,
        U=nothing,
        Q=nothing,
        T=nothing,
        F=nothing,
        L=nothing,
        ) = new(
            eta,
            D === nothing ? (u,N) -> u[id.D(N)] : D,
            C === nothing ? (u,N) -> u[id.C(N)] : C,
            R === nothing ? (u,N) -> u[id.R(N)] : R,
            H === nothing ? (u,N) -> u[id.H(N)] : H,
            U === nothing ? (u,N) -> u[id.U(N)] : U,
            Q === nothing ? (u,N) -> u[id.Q(N)] : Q,
            T === nothing ? (u,N) -> u[id.T(N)] : T,
            F === nothing ? (u,N) -> u[id.F(N)] : F,
            L === nothing ? (u,N) -> u[id.L(N)] : L,
        )

end

module Elevation
using ..Index:elevation_indexes

struct ElevationStruct
    point
    lowest
    highest
    average
    potential_energy
end

module Direct
    using ..Elevation: ElevationStruct, elevation_indexes

    function point(u,N,m)
        return u[elevation_indexes(N)[m+begin]]
    end

    function lowest(u,N)
        return u[elevation_indexes(N)[end]]
    end

    function highest(u,N)
        return u[elevation_indexes(N)[begin]]
    end

    function average(u,N)
        correction = (highest(u,N) + lowest(u,N)) / 2
        return (sum(u[elevation_indexes(N)]) - correction) / N
    end

    function potential_energy(u,N)
        relative_elevation = u[elevation_indexes(N)] .- average(u,N)

        return (relative_elevation[begin]^2 + relative_elevation[end]^2 + 2 * sum(relative_elevation[begin+1:end-1] .^ 2)) / 4N
    end

    DIRECT_ELEVATION = ElevationStruct(
        point,
        lowest,
        highest,
        average,
        potential_energy
    )
end

end

export C_INDEX, D_INDEX, H_INDEX, Q_INDEX, R_INDEX, U_INDEX

export elevation_indexes,stream_indexes

end
