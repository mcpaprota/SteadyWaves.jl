
module Surface
using ..Index: IndexStruct
using ..Params

struct SurfaceStruct
    point
    keypoints
    min
    max
    avg
    e_p
    a
    z

    SurfaceStruct(point,keypoints,min,max,avg,e_p,a,z) = new(point,keypoints,min,max,avg,e_p,a,z)

    SurfaceStruct(u,compiler::SurfaceStruct,wave_compiler) = new(
        compiler.point(wave_compiler,u),
        compiler.keypoints(wave_compiler,u),
        compiler.min(u),
        compiler.max(u),
        compiler.avg(u),
        compiler.e_p(wave_compiler,u),
        compiler.a(u),
        compiler.z(u)
    )

    SurfaceStruct(
        default::SurfaceStruct;
        point=nothing,
        keypoints=nothing,
        min=nothing,
        max=nothing,
        avg=nothing,
        e_p=nothing,
        a=nothing,
        z=nothing,
    ) = new(
            something(point,     default.point),
            something(keypoints, default.keypoints),
            something(min,       default.min),
            something(max,       default.max),
            something(avg,       default.avg),
            something(e_p,       default.e_p),
            something(a,         default.a),
            something(z,         default.z),
        )
    
    SurfaceStruct(idx::IndexStruct,eta_type::ElevationType) = begin
        if eta_type == Params.DIRECT_ELEVATION
            return direct_elevation_struct(idx)
        elseif eta_type == Params.FOURIER_ELEVATION
            return fourier_elevation_struct(idx)
        end
    end

    (compiler::SurfaceStruct)(w_c,u::Vector) = return SurfaceStruct(u,compiler,w_c)
end

struct EtaSupportStruct
    z
    x
    dz_dx_1
    dz_dx_2


    EtaSupportStruct(z,x,dz_dx_1,dz_dx_2) = new(z,x,dz_dx_1,dz_dx_2)


    (e_c::EtaSupportStruct)(w_c,u) = return EtaSupportStruct(
        e_c.z(w_c,u),
        e_c.x(w_c,u),
        e_c.dz_dx_1(w_c,u),
        e_c.dz_dx_2(w_c,u),
    )
    (e_c::EtaSupportStruct)(m) = return e_c.z(m)
end

function direct_point_der_1(point,m,idx::IndexStruct)
    dkx = pi / idx.N
    return (point(m-1) - point(m+1))/ 2dkx
end

function direct_point_der_2(point,m,idx::IndexStruct)
    dkx = pi / idx.N
    return (point(m-1) - 2point(m) + point(m+1)) / dkx^2
end

function direct_elevation_struct(idx)
    point = EtaSupportStruct(
        (w_c,u) -> m -> direct_point(u,m,idx),
        (w_c,u) -> m -> point_x(m,idx),
        (w_c,u) -> m -> direct_point_der_1(j -> direct_point(u,j,idx),m,idx),
        (w_c,u) -> m -> direct_point_der_2(j -> direct_point(u,j,idx),m,idx),
    )


    return SurfaceStruct(
        point,
        (w_c,u) -> direct_keypoints(u,idx),
        (u) -> direct_min(u,idx),
        (u) -> direct_max(u,idx),
        (u) -> direct_avg(u,idx),
        (w_c,u) -> direct_potential_energy(u,idx,w_c),
        u -> nothing,
        u -> nothing,
    )
end

function fourier_elevation_struct(idx)
    point = EtaSupportStruct(
        (w_c,u) -> m -> fourier_point(w_c,u,m,idx),
        (w_c,u) -> m -> point_x(m,idx),
        (w_c,u) -> m -> fourier_point_der_1(w_c,u,m,idx),
        (w_c,u) -> m -> fourier_point_der_2(w_c,u,m,idx),
    )

    return SurfaceStruct(
        point,
        (w_c,u) -> nothing,
        u -> fourier_min(u,idx),
        u -> fourier_max(u,idx),
        u -> fourier_avg(u,idx),
        (w_c,u) -> fourier_potential_energy(u,idx,w_c),
        u -> fourier_amplitudes(u,idx),
        u -> fourier_z_support_struct(u[idx.eta]),
    )
end

function point_x(m,idx::IndexStruct)
    return m * pi/idx.N
end

function keypoint_point(keypoint,m,idx::IndexStruct)
    println(m)
    m = idx.N - abs(idx.N - m%2idx.N)
    return keypoint[begin+m]
end

function direct_point(u,m,idx::IndexStruct)
    # triangular index function 
    P = 2idx.N
    m = (P + m%P)%P #ensure that m is in <0,P-1> range even if m < 0
    m = idx.N - abs(idx.N - m)
    return u[idx.eta[begin+m]]
end

function direct_keypoints(u,idx::IndexStruct)
    return u[idx.eta]
end

function direct_min(u,idx::IndexStruct)
    return u[idx.eta[end]]
end

function direct_max(u,idx::IndexStruct)
    return u[idx.eta[begin]]
end

function direct_avg(u,idx::IndexStruct)
    return (2 * sum(u[idx.eta]) - direct_min(u,idx) - direct_max(u,idx)) / 2 / idx.N 
end

function direct_potential_energy(u,idx::IndexStruct,w_c)
    relative_eta = u[idx.eta] .- w_c.D(w_c,u)
    
    return (relative_eta[1]^2 + relative_eta[end]^2 + 2 * sum(relative_eta[2:end-1] .^ 2)) / 4idx.N
end

function fourier_amplitudes(u,idx::IndexStruct)
    return u[idx.eta]
end

function fourier_z(a,kx)
    return sum(map(j -> a[begin+j] * cos(j * kx), 0:lastindex(a)-1))
end

function fourier_dz_dx_1(a,kx)
    return sum(map(j -> -j*a[begin+j] * sin(j * kx), 0:lastindex(a)-1))
end

function fourier_dz_dx_2(a,kx)
    return sum(map(j -> -(j^2)*a[begin+j] * cos(j * kx), 0:lastindex(a)-1))
end

function fourier_z(u,kx,idx)
    return fourier_z(u[idx.eta],kx)
end

function fourier_point(w_c,u,m::Int,idx::IndexStruct)
    return fourier_z(
        u[idx.eta],
        w_c.eta.point(w_c,u).x(m),
    )
end

function fourier_point_der_1(w_c,u,m,idx)
    return fourier_dz_dx_1(
        u[idx.eta],
        w_c.eta.point(w_c,u).x(m),
    )
end

function fourier_point_der_2(w_c,u,m,idx)
    return fourier_dz_dx_2(
        u[idx.eta],
        w_c.eta.point(w_c,u).x(m),
    )
end

function fourier_max(u,idx::IndexStruct)
    return fourier_z(u,0,idx)
end

function fourier_min(u,idx::IndexStruct)
    return fourier_z(u,pi,idx)
end

function fourier_avg(u,idx::IndexStruct)
    return u[idx.eta[begin]]
end

function fourier_potential_energy(u,idx,w_c)
    relative_a_0 = u[idx.eta[begin]] - w_c.D(w_c,u)

    return pi*relative_a_0^2  + sum(map(a -> pi * a^2, u[idx.eta[begin+1:end]]))
end

function fourier_from_points(point_z,point_x,N)
    a = zeros(N+1)
    
    a[begin]= (2 * sum(point_z.(0:N))-point_z(0)- point_z(N))/(2 * N)
    
    for j in 1:N
        coef = point_z.(0:N) .* cos.(j*point_x.(0:N))

        a[j+begin] = (2*sum(coef) - coef[begin] - coef[end])/N
    end
    return a
end

function struct_with_derived_values(eta::SurfaceStruct,idx::IndexStruct,type::ElevationType)
    if type == Params.DIRECT_ELEVATION
        return struct_with_fourier(eta,idx)
    elseif type == Params.FOURIER_ELEVATION
        return struct_with_keypoints(eta,idx)
    end
end

function fourier_z_support_struct(a)
    return EtaSupportStruct(
        kx -> fourier_z(a,kx),
        kx -> kx,
        kx -> fourier_dz_dx_1(a,kx),
        kx -> fourier_dz_dx_2(a,kx),
    )
    
end

function struct_with_fourier(eta::SurfaceStruct,idx::IndexStruct)
    a = fourier_from_points(eta.point,eta.point.x,idx.N)

    return SurfaceStruct(eta;
        a = a,
        z = fourier_z_support_struct(a),
    )
end

function struct_with_keypoints(eta::SurfaceStruct,idx::IndexStruct)
    return SurfaceStruct(eta;
        keypoints = eta.point.(0:idx.N)
    )
end

end