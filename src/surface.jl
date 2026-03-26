
module Surface
using ..Index: IndexStruct

struct SurfaceStruct
    point
    min
    max
    avg
    e_p

    SurfaceStruct(point,min,max,avg,e_p) = new(point,min,max,avg,e_p)

    SurfaceStruct(u,compiler::SurfaceStruct,wave_compiler) = new(
        m -> compiler.point(u,m),
        compiler.min(u),
        compiler.max(u),
        compiler.avg(u),
        compiler.e_p(wave_compiler,u)
    )

    SurfaceStruct(idx::IndexStruct) = new(
        (u,m) -> direct_point(u,m,idx),
        (u) -> direct_min(u,idx),
        (u) -> direct_max(u,idx),
        (u) -> direct_avg(u,idx),
        (w_c,u) -> potential_energy(u,idx,w_c)
    )

    (compiler::SurfaceStruct)(w_c,u::Vector) = return SurfaceStruct(u,compiler,w_c)

end

function direct_point(u,m,idx::IndexStruct)
    return u[idx.eta[begin+m]]
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

end