
module Surface
using ..Index: IndexStruct

struct SurfaceStruct
    point
    min
    max
    avg

    SurfaceStruct(point,min,max,avg) = new(point,min,max,avg)

    SurfaceStruct(u,compiler::SurfaceStruct) = new(
        m -> compiler.point(u,m),
        compiler.min(u),
        compiler.max(u),
        compiler.avg(u),
    )

    SurfaceStruct(idx::IndexStruct) = new(
        (u,m) -> direct_point(u,m,idx),
        (u) -> direct_min(u,idx),
        (u) -> direct_max(u,idx),
        (u) -> direct_avg(u,idx),
    )

    (compiler::SurfaceStruct)(u::Vector) = return SurfaceStruct(u,compiler)

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

end