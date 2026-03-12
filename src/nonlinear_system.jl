module NonlinearSystem

using ..Output
using ..Output: wave_power, wave_period, surface_stream_eigenfunction, wave_length
using ..Params

function mean_depth_condition(u,N)
    elevation = u[elevation_indexes(N)]
    return ((elevation[1] + elevation[end]) / 2 + sum(elevation[2:end-1])) / N - u[2N+D_INDEX]
end

function kinematic_surface_condition(u,N,m)
    stream_coeffs = @view u[stream_indexes(N)]
    depth = u[2N+D_INDEX]
    kx = m/N * π
    kz = u[m+1] # surface elevation
    Σ₁ = sum([surface_stream_eigenfunction(sinh, cos, stream_coeffs, depth, kx, kz, j) for j in 1:N])
    return Σ₁ - u[2N+U_INDEX] * (kz - depth) - u[2N+Q_INDEX]
end

function dynamic_surface_condition(u,N,m)
    stream_coeffs = @view u[stream_indexes(N)]
    depth = u[2N + D_INDEX] 
    kx = m/N * π
    kz_surface = u[m+1]
    Σ₂ = sum([j*surface_stream_eigenfunction(cosh,cos,stream_coeffs,depth,kx,kz_surface,j) for j in 1:N])
    Σ₃ = sum([j*surface_stream_eigenfunction(sinh,sin,stream_coeffs,depth,kx,kz_surface,j) for j in 1:N])
    return (-u[2N+U_INDEX] + Σ₂)^2 / 2 + Σ₃^2 / 2 + kz_surface - depth - u[2N+R_INDEX]
end

function height_condition(u,N,p)
    return u[1] - u[N+1] - u[2N+3] * p
end

function height_condition(u,N)
    @assert length(u) >= 2N+7 "height_condition require u[2N+7] or parameter p { H / d * m / M } to exist"
    return u[1] - u[N+1] - u[2N+7]
end

function power_condition(u,N,p)
    return wave_power(u, N) - p * √u[2N+3]^5
end

function euler_condition(u,N)
    return u[2N+U_INDEX] - u[2N+C_INDEX]
end

function stokes_condition(u,N)
    return euler_condition(u,N) - u[2N+Q_INDEX] / u[2N+D_INDEX]
end

function length_condition(u,N,p)
    return u[2N+D_INDEX] - 2π / p
end

function period_condition(u,N,p)
    return u[2N+C_INDEX] * p * √u[2N+D_INDEX] - 2π   
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

function nonlinear_system_base!(du, u, N, conditions)

    index = 1

    for con_struct in conditions
        if con_struct.is_singular
            du[index] = con_struct.condition(u,N)
            index += 1
        else
            for m in con_struct.range
                du[index] = con_struct.condition(u,N,m)
                index += 1
            end
        end

    end

    return nothing
end

end