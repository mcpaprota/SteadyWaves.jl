function stream_eigenfunction(hiperbolic,trigonometric,u,N,m,j)
    return u[N+1+j] * hiperbolic(j * u[m+1]) / cosh(j * u[2N+3]) * trigonometric(j * m * π / N)
end

function mean_depth(u,N)
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


include("params.jl") 
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