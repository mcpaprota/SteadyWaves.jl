
module Surface
using ..Index: IndexStruct

struct SurfaceStruct
    point
    point_der_1
    point_der_2
    min
    max
    avg
    e_p

    SurfaceStruct(point,min,max,avg,e_p) = new(point,min,max,avg,e_p)

    SurfaceStruct(u,compiler::SurfaceStruct,wave_compiler) = new(
        m -> compiler.point(u,m),
        m -> compiler.point_der_1(wave_compiler,u,m),
        m -> compiler.point_der_2(wave_compiler,u,m),
        compiler.min(u),
        compiler.max(u),
        compiler.avg(u),
        compiler.e_p(wave_compiler,u)
    )

    SurfaceStruct(idx::IndexStruct) = new(
        (u,m) -> direct_point(u,m,idx),
        (w_c,u,m) -> direct_point_der_1(m -> w_c.eta.point(u,m),m,idx),
        (w_c,u,m) -> direct_point_der_2(m -> w_c.eta.point(u,m),m,idx),
        (u) -> direct_min(u,idx),
        (u) -> direct_max(u,idx),
        (u) -> direct_avg(u,idx),
        (w_c,u) -> potential_energy(u,idx,w_c)
    )

    (compiler::SurfaceStruct)(w_c,u::Vector) = return SurfaceStruct(u,compiler,w_c)

end

function direct_point(u,m,idx::IndexStruct)
    P = 2idx.N
    m = (P + m%P)%P #ensure that m is in <0,P-1> range even if m < 0
    m = idx.N - abs(idx.N - m)
    return u[idx.eta[begin+m]]
end

function direct_point_der_1(point,m,idx::IndexStruct)
    dkx = pi / idx.N
    return (point(m-1) - point(m+1))/ 2dkx
end

function direct_point_der_2(point,m,idx::IndexStruct)
    dkx = pi / idx.N
    return (point(m-1) - 2point(m) + point(m+1)) / dkx^2
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

function fourier_point(u,m,idx::IndexStruct)
    return fourier_z(u,pi * m/idx.N, idx)
end

function fourier_z(u,kx::Float64,idx::IndexStruct)
    return u[idx.D] + sum(map(j -> a_f[j] * cos(j * kx), idx.eta))
end

function fourier_z(a_f,kx::Float64,kd::Float64)
    return kd + sum(map(j -> a_f[j] * cos(j * kx), 1:lastindex(a_f)))
end

end