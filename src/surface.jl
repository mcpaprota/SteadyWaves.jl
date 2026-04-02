
module Surface
using ..Index: IndexStruct

struct SurfaceStruct
    point
    point_x
    keypoints
    min
    max
    avg
    e_p
    a
    z

    SurfaceStruct(point,point_x,min,max,avg,e_p,a,z) = new(point,point_x,min,max,avg,e_p,a,z)

    SurfaceStruct(u,compiler::SurfaceStruct,wave_compiler) = new(
        m -> compiler.point(u,m),
        m -> compiler.point_x(u,m),
        compiler.keypoints(u),
        compiler.min(u),
        compiler.max(u),
        compiler.avg(u),
        compiler.e_p(wave_compiler,u),
        compiler.a(u),
        kx -> compiler.z(u,kx)
    )

    SurfaceStruct(
        default::SurfaceStruct;
        point=nothing,
        point_x= nothing,
        keypoints=nothing,
        min=nothing,
        max=nothing,
        avg=nothing,
        e_p=nothing,
        a=nothing,
        z=nothing,
    ) = new(
            something(point,     default.point),
            something(point_x,   default.point_x),
            something(keypoints, default.keypoints),
            something(min,       default.min),
            something(max,       default.max),
            something(avg,       default.avg),
            something(e_p,       default.e_p),
            something(a,         default.a),
            something(z,         default.z),
        )
    
    SurfaceStruct(idx::IndexStruct) = new(
        (u,m) -> direct_point(u,m,idx),
        (u,m) -> direct_point_x(u,m,idx),
        (u) -> direct_keypoints(u,idx),
        (u) -> direct_min(u,idx),
        (u) -> direct_max(u,idx),
        (u) -> direct_avg(u,idx),
        (w_c,u) -> potential_energy(u,idx,w_c),
        u -> nothing,
        u -> nothing,
    )



    (compiler::SurfaceStruct)(w_c,u::Vector) = return SurfaceStruct(u,compiler,w_c)

end




function direct_point(u,m,idx::IndexStruct)
    # triangular index function 
    m = idx.N - abs(idx.N - m%2idx.N)
    return u[idx.eta[begin+m]]
end

function direct_point_x(u,m,idx::IndexStruct)
    return m * pi/idx.N
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

function potential_energy(u,idx::IndexStruct,w_c)
    relative_eta = u[idx.eta] .- w_c.D(w_c,u)
    
    return (relative_eta[1]^2 + relative_eta[end]^2 + 2 * sum(relative_eta[2:end-1] .^ 2)) / 4idx.N
end

function fourier_z(a,kx)
    return sum(map(j -> a[begin+j] * cos(j * kx), 0:lastindex(a)-1))
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

function struct_with_fourier(eta::SurfaceStruct,idx::IndexStruct)
    a = fourier_from_points(eta.point,eta.point_x,idx.N)

    return SurfaceStruct(eta;
        a = a,
        z = kx -> fourier_z(a,kx),
    )
end

end