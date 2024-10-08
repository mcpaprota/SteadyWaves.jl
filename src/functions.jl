"""
    fourier_approx(d, H, L; cc=1, N=10, M=1)

Approximate wave solution `u` using Fourier Approximation Method.

...
# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `L`: wavelength (m)
- `cc`: current criterion 1 - Stokes; 2 - euler
- `N`: number of solution eigenvalues
- `M`: number of height

# Output (nondimensional)
- `u[1:N+1]`: free surface elevation *kη*
- `u[N+2:2N+1]`: stream function coefficients *B*
- `u[2N+2]`: wave celerity *c√(k/g)*
- `u[2N+3]`: mean water depth *kη̄*
- `u[2N+4]`: volume flux due to waves *q√(k³/g)*
- `u[2N+5]`: Bernoulli constant *rk/g*
- `u[2N+6]`: mean flow velocity *Ū√(k/g)*
"""
function fourier_approx(d, H, L; cc=1, N=10, M=1)
    u = init_conditions(d, H, L, N, M)
    for m in 1:M
        params = [L / d, H / d * m / M, cc]
        problem = NonlinearProblem(f, u, params)
        solution = solve(problem, RobustMultiNewton())
        u[:] = solution.u
    end
    return u
end

function init_conditions(d, H, L, N, M)
    k = 2π / L # wave number (rad/m)
    u0 = zeros(2N + 6)
    u0[1:N+1] = @. k * d + 1 / 2 * k * H / M * cos((0:N) * π / N) # kη
    u0[N+2:2N+1] = [k * H / M / 2 / √tanh(k * d); zeros(N - 1)] # B
    u0[2N+2] = √tanh(k * d) # c√(k/g)
    u0[2N+3] = k * d # kη̄
    u0[2N+4] = 0 # q√(k³/g)
    u0[2N+5] = tanh(k * d) / 2 # rk/g
    u0[2N+6] = √tanh(k * d) # Ū√(k/g)
    return u0
end

function f(du, u, p)
    N = (length(u) - 6) ÷ 2
    for m in 0:N
        Σ₁ = sum([u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N) for j in 1:N])
        du[m+1] = Σ₁ - u[2N+6] * (u[m+1] - u[2N+3]) - u[2N+4]
        Σ₂ = sum([j * u[N+1+j] * cosh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N) for j in 1:N])
        Σ₃ = sum([j * u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * sin(j * m * π / N) for j in 1:N])
        du[N+1+m+1] = (-u[2N+6] + Σ₂)^2 / 2 + Σ₃^2 / 2 + u[m+1] - u[2N+3] - u[2N+5]
    end
    du[2N+3] = ((u[1] + u[N+1]) / 2 + sum(u[2:N])) / N - u[2N+3]
    du[2N+4] = u[1] - u[N+1] - u[2N+3] * p[2]
    du[2N+5] = u[2N+3] - 2π / p[1]
    p[3] == 1 ? du[2N+6] = u[2N+6] - u[2N+2] - u[2N+4] / u[2N+3] : du[2N+6] = u[2N+6] - u[2N+2]
    return nothing
end
