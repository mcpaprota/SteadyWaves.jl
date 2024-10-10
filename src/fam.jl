# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)

"""
    fourier_approx(d, H, P; pc=1, cc=1, N=10, M=1, g=9.81)

Approximate solution `u` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using Fourier Approximation Method.

...
# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=1` - length (default), `pc=2` - period
- `cc`: current criterion; `cc=1` - Stokes (default), `cc=2` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
- `M`: number of height steps, defaults to `M=1`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`

# Output
- `u[1:N+1]`: free surface elevation *kη*
- `u[N+2:2N+1]`: stream function coefficients *B*
- `u[2N+2]`: wave celerity *c√(k/g)*
- `u[2N+3]`: mean water depth *kη̄*
- `u[2N+4]`: volume flux due to waves *q√(k³/g)*
- `u[2N+5]`: Bernoulli constant *rk/g*
- `u[2N+6]`: mean flow velocity *Ū√(k/g)*
"""
function fourier_approx(d, H, P; pc=1, cc=1, N=10, M=1, g=9.81)
    u = init_conditions(d, H, P, pc, N, M)
    for m in 1:M
        if pc==1
            params = [P / d, H / d * m / M, pc, cc]
        else
            params = [P * sqrt(g / d), H / d * m / M, pc, cc]
        end
        problem = NonlinearProblem(f_0, u, params)
        solution = solve(problem, RobustMultiNewton())
        u[:] = solution.u
    end
    return u
end

"""
    init_conditions(d, H, P, pc, N, M)

Calculate initial conditions `u0` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using linear wave theory.

...
# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=1` - length, `pc=2` - period
- `N`: number of solution eigenvalues
- `M`: number of height steps

# Output
- `u[1:N+1]`: free surface elevation *kη*
- `u[N+2:2N+1]`: stream function coefficients *B*
- `u[2N+2]`: wave celerity *c√(k/g)*
- `u[2N+3]`: mean water depth *kη̄*
- `u[2N+4]`: volume flux due to waves *q√(k³/g)*
- `u[2N+5]`: Bernoulli constant *rk/g*
- `u[2N+6]`: mean flow velocity *Ū√(k/g)*
"""
function init_conditions(d, H, P, pc, N, M)
    pc == 1 ? k = 2π / P : k = dispersion_relation(d, 2π / P) # wave number (rad/s)
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

"""
    f_0(du, u, p)

Define nonlinear system `f_0(u) = 0` with parameters `p`.

"""
function f_0(du, u, p)
    N = (length(u) - 6) ÷ 2
    for m in 0:N
        Σ₁ = sum([u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N)
                  for j in 1:N])
        du[m+1] = Σ₁ - u[2N+6] * (u[m+1] - u[2N+3]) - u[2N+4]
        Σ₂ = sum([j * u[N+1+j] * cosh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N)
                  for j in 1:N])
        Σ₃ = sum([j * u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * sin(j * m * π / N)
                  for j in 1:N])
        du[N+1+m+1] = (-u[2N+6] + Σ₂)^2 / 2 + Σ₃^2 / 2 + u[m+1] - u[2N+3] - u[2N+5]
    end
    du[2N+3] = ((u[1] + u[N+1]) / 2 + sum(u[2:N])) / N - u[2N+3]
    du[2N+4] = u[1] - u[N+1] - u[2N+3] * p[2]
    if p[3] == 1
        du[2N+5] = u[2N+3] - 2π / p[1]
    else
        du[2N+5] = u[2N+2] * p[1] * √u[2N+3] - 2π
    end
    if p[4] == 1
        du[2N+6] = u[2N+6] - u[2N+2] - u[2N+4] / u[2N+3]
    else
        du[2N+6] = u[2N+6] - u[2N+2]
    end
    return nothing
end

"""
    dispertion_relation(d, ω, g=9.81, ϵ=10^-12)

Calculate wavenumber `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""
function dispersion_relation(d, ω, g=9.81, ϵ=10^-12)
    k = k₀ = ω^2 / g # initial guess
    while max(abs(k * tanh(k * d) - k₀)) > ϵ
        k = k₀ / tanh(k * d)
    end
    return k
end
