module Velocity

using ..Params
    
struct VelocityStruct
    x
    z
    psi

    VelocityStruct(x,z,psi) = new(x,z,psi)

    (v::VelocityStruct)(w_c,u) = return VelocityStruct(
        v.x(w_c,u),
        v.z(w_c,u),
        v.psi(w_c,u)
    )
end

function shallow_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, shallow_stream_eigenfunction(cosh,cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, shallow_stream_eigenfunction(sinh,sin)),
        (w_c,u) -> stream(                      u, idx, w_c, shallow_stream_eigenfunction(sinh,cos)),
    )
end

function deep_water_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, deep_water_stream_eigenfunction(cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, deep_water_stream_eigenfunction(sin)),
        (w_c,u) -> stream(                      u, idx, w_c, deep_water_stream_eigenfunction(cos)),
    )
end

function stable_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, stable_stream_eigenfunction(cosh,sinh,cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, stable_stream_eigenfunction(sinh,cosh,sin)),
        (w_c,u) -> stream(                      u, idx, w_c, stable_stream_eigenfunction(sinh,cosh,cos)),
    )
end

function velocity_struct_factory(idx,config)
    if config.deep_water == Params.DEEP_WATER
        return deep_water_velocity_struct(idx)
    elseif config.deep_water == Params.SHALLOW_WATER
        return shallow_velocity_struct(idx)
    else
        return stable_velocity_struct(idx)
    end
end


function shallow_stream_eigenfunction(hiperbolic,trigonometric)
    return (B,kd,kx,kz,j) -> B * hiperbolic(j * kz) / cosh(j * kd) * trigonometric(j * kx)
end

function stable_stream_eigenfunction(h_single,h_multiplied,trigonometric)
    return (B,kd,kx,kz,j) -> B * (h_single(j * (kz-kd)) + tanh(j * kd)*h_multiplied(j * (kz-kd)) )* trigonometric(j * kx)
end

function stable_stream_eigenfunction(multi,trigonometric)
    return (B,kd,kx,kz,j) -> B * (exp(kz-kd) + multi * exp(-kz-kd)) / (1 + exp(-kd)) * trigonometric(j * kx)
end


function deep_water_stream_eigenfunction(trigonometric)
    return (B,kd,kx,kz,j) -> B * exp(j*(kz-kd))*trigonometric(j*kx)
end


function stream_horizontal_velocity(u, idx, w_c,stream_eigenfunction)
    kd = w_c.D(w_c,u)
    return (kx,kz) -> -w_c.U(w_c,u) + sum([j*stream_eigenfunction(u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream_vertical_velocity(u, idx, w_c,stream_eigenfunction)
    kd = w_c.D(w_c,u)
    return (kx,kz) -> sum([j*stream_eigenfunction(u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream(u, idx, w_c, stream_eigenfunction)
    kd = w_c.D(w_c,u)
    return (kx,kz) -> sum([stream_eigenfunction(u[idx.psi[j]], kd, kx, kz, j) for j in 1:idx.N])
end

end