module Velocity
    
struct VelocityStruct
    x
    z
    psi

    VelocityStruct(x,z,psi) = new(x,z,psi)

    VelocityStruct(idx) = new(
        (w_c,u) -> (kx,kz) -> stream_horizontal_velocity(u,idx,w_c,kx,kz),
        (w_c,u) -> (kx,kz) -> stream_vertical_velocity(u,idx,w_c,kx,kz),
        (w_c,u) -> (kx,kz) -> stream(u,idx,w_c,kx,kz),
    )

    (v::VelocityStruct)(w_c,u) = return VelocityStruct(
        v.x(w_c,u),
        v.z(w_c,u),
        v.psi(w_c,u)
    )
end

function stream_eigenfunction(hiperbolic,trigonometric,B,kd,kx,kz,j)
    return B * hiperbolic(j * kz) / cosh(j * kd) * trigonometric(j * kx)
end

function stream_horizontal_velocity(u, idx, w_c, kx, kz)
    kd = w_c.D(w_c,u)
    return -w_c.U(w_c,u) + sum([j*stream_eigenfunction(cosh, cos, u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream_vertical_velocity(u, idx, w_c, kx, kz)
    kd = w_c.D(w_c,u)
    return sum([j*stream_eigenfunction(sinh, sin, u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream(u, idx, w_c, kx, kz)
    kd = w_c.D(w_c,u)
    return sum([stream_eigenfunction(sinh, cos, u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

end