# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)
include("params.jl")
include("nonlinear_system.jl")
"""
    fourier_approx(d, H, P; pc=1, cc=1, N=10, M=1, g=9.81)

Approximate solution `u` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using Fourier Approximation Method.

# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=1`, `pc=PC_LENGTH` - length (default), `pc=2`, `pc=PC_PERIOD` - period
- `cc`: current criterion; `cc=1`, `cc=CC_STOKES` - Stokes (default), `cc=2`, `cc=CC_EULER` - Euler
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
- `u[2N+7]`: wave height *kH*
"""
function fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=9.81)
    u = init_conditions(d, H, P, Int(pc), N, M)
    for m in 1:M
        if Int(pc) == Int(PC_LENGTH)
            params = [P / d, H / d * m / M, Int(pc), Int(cc)]
        elseif Int(pc) == Int(PC_PERIOD)
            params = [P * sqrt(g / d), H / d * m / M, Int(pc), Int(cc)]
        else # prevents parsing unknown parameters
            throw(error("Unknown parameter criterion \"$pc\""))
        end

        problem = NonlinearProblem(nonlinear_system_steady, u, params)
        solution = solve(problem, RobustMultiNewton())
        u[:] = solution.u
    end
    push!(u, u[2N+3] / d * H)
    return u
end

"""
    init_conditions(d, H, P, pc, N, M)

Calculate initial conditions `u0` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using linear wave theory.

# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=1`, `pc=PC_LENGTH` - length (default), `pc=2`, `pc=PC_PERIOD` - period
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
    Int(pc) == Int(PC_LENGTH) ? k = 2π / P : k = wave_number(d, 2π / P) # wave number (rad/s)
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
    nonlinear_system_steady(du, u, p)

Define nonlinear system for steady waves `f(u) = 0` with parameters `p`.

"""
function nonlinear_system_steady(du, u, p)
    N = (length(u) - 6) ÷ 2
    
    for m in 0:N
        du[m+1] = kinematic_surface_condition(u, N, m)
        du[N+1+m+1] = dynamic_surface_condition(u, N, m)
    end

    du[2N+3] = mean_depth(u,N)
    
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
    wave_number(d, ω, g=9.81, ϵ=10^-12)

Calculate wave_number `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""
function wave_number(d, ω, g=9.81, ϵ=10^-12)
    k = k₀ = ω^2 / g # initial guess
    while max(abs(k * tanh(k * d) - k₀)) > ϵ
        k = k₀ / tanh(k * d)
    end
    return k
end
