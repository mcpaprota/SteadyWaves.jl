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

end