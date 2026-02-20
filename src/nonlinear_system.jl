include("output.jl")
include("params.jl") 

function stream_eigenfunction(hiperbolic,trigonometric,u,N,m,j)
    return u[N+1+j] * hiperbolic(j * u[m+1]) / cosh(j * u[2N+3]) * trigonometric(j * m * π / N)
end

function mean_depth_condition(u,N)
    return ((u[1] + u[N+1]) / 2 + sum(u[2:N])) / N - u[2N+3]
end

function kinematic_surface_condition(u,N,m)
    Σ₁ = sum([stream_eigenfunction(sinh, cos, u, N, m, j) for j in 1:N])
    return Σ₁ - u[2N+6] * (u[m+1] - u[2N+3]) - u[2N+4]
end

function dynamic_surface_condition(u,N,m)
        Σ₂ = sum([j*stream_eigenfunction(cosh,cos,u,N,m,j) for j in 1:N])
        Σ₃ = sum([j*stream_eigenfunction(sinh,sin,u,N,m,j) for j in 1:N])
        return (-u[2N+6] + Σ₂)^2 / 2 + Σ₃^2 / 2 + u[m+1] - u[2N+3] - u[2N+5]
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
    return u[2N+6] - u[2N+2]
end

function stokes_condition(u,N)
    return euler_condition(u,N) - u[2N+4] / u[2N+3]
end

function length_condition(u,N,p)
    return u[2N+3] - 2π / p
end

function period_condition(u,N,p)
    return u[2N+2] * p * √u[2N+3] - 2π   
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

function nonlinear_system_base(du, u, N, h_condition, pc_equation, cc_equation, pwr_condition)
    for m in 0:N
        du[m+1] = kinematic_surface_condition(u, N, m)
        du[N+1+m+1] = dynamic_surface_condition(u, N, m)
    end

    du[2N+3] = mean_depth_condition(u,N)
    
    du[2N+4] = h_condition(u,N)

    du[2N+5] = pc_equation(u,N)

    du[2N+6] = cc_equation(u,N)

    if !isnothing(pwr_condition)
        du[2N+7] = pwr_condition(u,N)
    end
    return nothing
end
