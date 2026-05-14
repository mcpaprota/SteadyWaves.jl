module Index

function eta_indexes(N)
return 1:N+1
end
function psi_indexes(N)
return N+2:2N+1
end
const C_INDEX::Int = 2
const D_INDEX::Int = 3
const Q_INDEX::Int = 4
const R_INDEX::Int = 5
const U_INDEX::Int = 6
const H_INDEX::Int = 7

struct IndexStruct
    eta::UnitRange
    psi::UnitRange
    D::Int
    C::Int
    R::Int
    H::Int
    U::Int
    Q::Int
    N::Int
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

function dynamic_indexes(N;eta=false, psi=false, D=false, C=false, R=false, H=false, U=false, Q=false)
    offset = 1

    eta, offset = dynamic_range(eta, offset)
    psi, offset = dynamic_range(psi, offset)
    D, offset   = dynamic_position(D, offset)
    C, offset   = dynamic_position(C, offset)
    R, offset   = dynamic_position(R, offset)
    H, offset   = dynamic_position(H, offset)
    U, offset   = dynamic_position(U, offset)
    Q, offset   = dynamic_position(Q, offset)

    return IndexStruct(eta, psi, D, C, R, H, U, Q,N)
end

function dynamic_range(values,offset)
    if values != false
        return (offset - values[begin]) .+ values, values[end] - values[begin] + offset + 1
    else
        return 0:0, offset
    end
end

function dynamic_position(value,offset)
    if value
        return offset, offset+1
    else
        return 0, offset
    end
end

function max_index(index::IndexStruct)
    return max(
        index.eta[end],
        index.psi[end],
        index.D,
        index.C,
        index.R,
        index.H,
        index.U,
        index.Q,
    )
end

end