module NonlinearSystem

using ..Index
using ..Output
using ..Output: wave_power, wave_period, surface_stream_eigenfunction, dimensionless_pressure
using ..Params

function mean_depth_condition(u,N,v)
    return v.eta.average(u,N) - v.D(u,N)
end

function kinematic_surface_condition(u,N,m,v)
    Σ₁ = sum([surface_stream_eigenfunction(sinh, cos, u, N, m, j) for j in 1:N])
    return Σ₁ - v.U(u,N) * (v.eta.point(u,N,m) - v.D(u,N)) - v.Q(u,N)
end

function dynamic_surface_condition(u,N,m,v)
    kx = m/N * pi
    kz = v.eta.point(u,N,m)

    return dimensionless_pressure(u, N, kx, kz)
end

function height_condition(u,N,H_d,v)
    return v.eta.highest(u,N) - v.eta.lowest(u,N) - v.D(u,N) * H_d
end

function height_condition(u,N,v)
    return v.eta.highest(u,N) - v.eta.lowest(u,N) - v.H(u,N)
end

function power_condition(u,N,v)
    return wave_power(u, N) - v.F(u,N)
end

function euler_condition(u,N,v)
    return v.U(u,N) - v.C(u,N)
end

function stokes_condition(u,N,v)
    return euler_condition(u,N,v) - v.U(u,N) / v.D(u,N)
end

function length_condition(u,N,v)
    return v.L(u,N) - 2π
end

function period_condition(u,N,v)
    return v.C(u,N) * v.T(u,N)  - 2π
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

function nonlinear_system_base!(du, u, N, conditions,v)

    index = 1

    for con_struct in conditions
        if con_struct.is_singular
            du[index] = con_struct.condition(u,N,v)
            index += 1
        else
            for m in con_struct.range
                du[index] = con_struct.condition(u,N,m,v)
                index += 1
            end
        end

    end

    return nothing
end

end