module Velocity

using ..Params
using ..StructOperator
"""
Structure with velocity field properties:
- `x`: horizontal velocity
- `z`: vertical velocity
- `psi`: value of stream function
- `b`: dimensionless amplitudes of stream function
"""
struct VelocityStruct
    x
    z
    psi
    b

    VelocityStruct(x,z,psi,b) = new(x,z,psi,b)

    (v_c::VelocityStruct)(w_c,u) = return StructOperator.map(v_c, v->v(w_c,u))
end

function shallow_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, shallow_stream_eigenfunction(cosh,cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, shallow_stream_eigenfunction(sinh,sin)),
        (w_c,u) -> stream(                      u, idx, w_c, shallow_stream_eigenfunction(sinh,cos)),
        b(idx),
    )
end

function deep_water_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, deep_water_stream_eigenfunction(cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, deep_water_stream_eigenfunction(sin)),
        (w_c,u) -> stream(                      u, idx, w_c, deep_water_stream_eigenfunction(cos)),
        b(idx),
    )
end

function stable_velocity_struct(idx)
    return VelocityStruct(
        (w_c,u) -> stream_horizontal_velocity(  u, idx, w_c, stable_stream_eigenfunction(cosh,sinh,cos)),
        (w_c,u) -> stream_vertical_velocity(    u, idx, w_c, stable_stream_eigenfunction(sinh,cosh,sin)),
        (w_c,u) -> stream(                      u, idx, w_c, stable_stream_eigenfunction(sinh,cosh,cos)),
        b(idx),
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
    return (kx,kz) -> -w_c.U(w_c,u) + sum([j*stream_eigenfunction(u[idx.v[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream_vertical_velocity(u, idx, w_c,stream_eigenfunction)
    kd = w_c.D(w_c,u)
    return (kx,kz) -> sum([j*stream_eigenfunction(u[idx.v[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function stream(u, idx, w_c, stream_eigenfunction)
    kd = w_c.D(w_c,u)
    println(typeof(u[idx.v]))
    return (kx,kz) -> sum([stream_eigenfunction(u[idx.v[j]], kd, kx, kz, j) for j in 1:idx.N])
end

function b(idx)
    return (w_c,u) -> u[idx.v]
end

end