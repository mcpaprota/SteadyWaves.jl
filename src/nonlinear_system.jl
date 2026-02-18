include("output.jl")
include("params.jl") 

function stream_eigenfunction(hiperbolic,trigonometric,u,N,m,j)
    return u[N+1+j] * hiperbolic(j * u[m+1]) / cosh(j * u[2N+D_INDEX]) * trigonometric(j * m * π / N)
end

function mean_depth_condition(u,N)
    elevation = u[elevation_indexes(N)]
    return ((elevation[1] + elevation[end]) / 2 + sum(elevation[2:end-1])) / N - u[2N+D_INDEX]
end

function kinematic_surface_condition(u,N,m)
    Σ₁ = sum([stream_eigenfunction(sinh, cos, u, N, m, j) for j in 1:N])
    return Σ₁ - u[2N+U_INDEX] * (u[m+1] - u[2N+D_INDEX]) - u[2N+Q_INDEX]
end

function dynamic_surface_condition(u,N,m)
        Σ₂ = sum([j*stream_eigenfunction(cosh,cos,u,N,m,j) for j in 1:N])
        Σ₃ = sum([j*stream_eigenfunction(sinh,sin,u,N,m,j) for j in 1:N])
        return (-u[2N+U_INDEX] + Σ₂)^2 / 2 + Σ₃^2 / 2 + u[m+1] - u[2N+D_INDEX] - u[2N+R_INDEX]
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

function current_criterion(cc)
    if Int(cc) == Int(CC_STOKES)
        return stokes_condition
    elseif Int(cc) == Int(CC_EULER)
        return euler_condition
    else
        throw(error("Unknown current criterion $cc"))
    end
end

function parameter_criterion(pc)
    if Int(pc) == Int(PC_LENGTH)
        return length_condition
    elseif Int(pc) == Int(PC_PERIOD)
        return period_condition
    else
        throw(error("Unknown parameter criterion $pc"))
    end
end

function parameter_criterion_constant(pc, P, d, g)
    if Int(pc) == Int(PC_LENGTH)
        return  P / d
    elseif Int(pc) == Int(PC_PERIOD)
        return P * sqrt(g / d) 
    else
        throw(error("Unknown parameter criterion $pc"))
    end
end