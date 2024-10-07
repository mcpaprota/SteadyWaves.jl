"""
    fourier_approx(d, H, L; cc=1, N=10, M=1)

Approximate wave solution `u` using Fourier Approximation Method.

...
# Arguments
- `d::Real`: water depth (m)
- `H::Real`: wave height (m)
- `L::Real`: wavelength (m)
- `cc:Integer`: current criterion 1 - Stokes; 2 - euler
- `N::Integer`: number of solution eigenvalues
- `M::Integer`: number of height

# Output (nondimensional)
- `u[1:N+1]`: free surface elevation *kη*
- `u[N+2:2N+1]`: stream function coefficients *B*
- `u[2N+2]`: wave celerity *c√(k/g)*
- `u[2N+3]`: mean water depth *kη̄*
- `u[2N+4]`: volume flux due to waves *q√(k³/g)*
- `u[2N+5]`: Bernoulli constant *rk/g*
- `u[2N+6]`: mean flow velocity *Ū√(k/g)*
"""
function fourier_approx(d::Real, H::Real, L::Real; cc::Integer=1, N::Integer=10, M::Integer=1)
    u0 = init_conditions(d, H, L, cc, N, M)
    for m in 1:M
        params = [L / d, H / d * m / M, cc]
        problem = NonlinearProblem(f, u0, params)
        solution = solve(problem, RobustMultiNewton)
        u[:] = solution.u
    end
    return u
end

function init_conditions(d, H, L, cc, N, M)
    k = 2π / L # wave number (rad/m)
    u0 = zeros(2N + 6)
    u0[1:N+1] = @. k * d + 1 / 2 * k * H / M * cos((0:N) * π / N) # kη
    u0[N+2:2N+1] = [k * H / M / 2 / √tanh(k * d); zeros(N - 1)] # B
    u0[2N+2] = √tanh(k * d) # c√(k/g)
    u0[2N+3] = k * d # kη̄
    u0[2N+4] = 0 # q√(k³/g)
    u0[2N+5] = tanh(k * d) / 2 # rk/g
    u0[2N+6] = √tanh(k * d) # Ū√(k/g)
end
