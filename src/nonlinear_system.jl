module NonlinearSystem

using ..Index:IndexStruct
using ..Output
using ..Output: wave_power, wave_period, surface_stream_eigenfunction, dimensionless_pressure
using ..Params

function mean_depth_condition(u,N,idx::IndexStruct)
    eta = u[idx.eta]
    return ((eta[begin] + eta[end]) / 2 + sum(eta[2:end-1])) / N - u[idx.D]
end

function kinematic_surface_condition(u,N,idx::IndexStruct,m)
    psi = sum([surface_stream_eigenfunction(sinh, cos, u, N, m, j) for j in 1:N])
    return psi - u[idx.U] * (u[idx.eta[begin+m]] - u[idx.D]) - u[idx.Q]
end

function dynamic_surface_condition(u,N,idx::IndexStruct,m)
    kx = m/N * pi
    kz = u[idx.eta[m+begin]]

    return dimensionless_pressure(u, N, kx, kz)
end

function height_condition(u,N,idx::IndexStruct,p)
    return u[idx.eta[begin]] - u[idx.eta[end]] - u[idx.D] * p
end

function height_condition(u,N,idx::IndexStruct)
    @assert length(u) >= 2N+7 "height_condition require u[2N+7] or parameter p { H / d * m / M } to exist"
    return u[idx.eta[begin]] - u[idx.eta[end]] - u[idx.H]
end

function power_condition(u,N,idx::IndexStruct,p)
    return wave_power(u, N) - p * √u[idx.D]^5
end

function euler_condition(u,N,idx::IndexStruct)
    return u[idx.U] - u[idx.C]
end

function stokes_condition(u,N,idx::IndexStruct)
    return euler_condition(u,N,idx::IndexStruct) - u[idx.Q] / u[idx.D]
end

function length_condition(u,N,idx::IndexStruct,p)
    return u[idx.D] - 2π / p
end

function period_condition(u,N,idx::IndexStruct,p)
    return u[idx.C] * p * √u[idx.D] - 2π   
end

function current_condition_factory(cc)
    if Int(cc) == Int(CC_STOKES)
        return stokes_condition
    elseif Int(cc) == Int(CC_EULER)
        return euler_condition
    else
        throw(error("Unknown current criterion $cc"))
    end
end

function parameter_condition_factory(pc)
    if Int(pc) == Int(PC_LENGTH)
        return length_condition
    elseif Int(pc) == Int(PC_PERIOD)
        return period_condition
    else
        throw(error("Unknown parameter criterion $pc"))
    end
end

function parameter_condition_constant(pc, P, d, g)
    if Int(pc) == Int(PC_LENGTH)
        return  P / d
    elseif Int(pc) == Int(PC_PERIOD)
        return P * sqrt(g / d) 
    else
        throw(error("Unknown parameter criterion $pc"))
    end
end

struct ConditionStruct
    condition
    range
    is_singular

    ConditionStruct(condition) = new(condition,nothing,true)

    ConditionStruct(condition,range) = new(condition,range,false)

end

function nonlinear_system_base!(du, u, N, conditions,idx::IndexStruct)

    index = 1

    for con_struct in conditions
        if con_struct.is_singular
            du[index] = con_struct.condition(u,N,idx)
            index += 1
        else
            for m in con_struct.range
                du[index] = con_struct.condition(u,N,idx,m)
                index += 1
            end
        end

    end

    return nothing
end

end