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