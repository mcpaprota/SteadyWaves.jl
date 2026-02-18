# SPDX-License-Identifier: MIT

# Functions for shoaling calculations based on Fourier Approximation Method
include("params.jl")
include("nonlinear_system.jl")

"""
    shoaling_approx(d, H, L; cc=2, N=10, g=9.81)

Calculate shoaling coefficients `K` in range of depth values `d`
for wave of length `L` and height `H`.

# Arguments
- `d`: vector of decreasing water depths (m)
- `L`: initial wavelength (m) - corresponding to d[1]
- `H`: initial wave height (m) - corresponding to d[1]
- `cc`: current criterion; `cc=1`, `cc=CC_STOKES` - Stokes (default), `cc=2`, `cc=CC_EULER` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`

# Output
- `K`: vector of shoaling coefficient values
"""
function shoaling_approx(d, H, L; cc=CC_STOKES, N=10, g=9.81)
    k = 2π / L # initial wave number (rad/m
    ω = √(g * k * tanh(k * d[1])) # initial angular wave frequency (rad/s)

    K = zero(float(d))
    K[1] = 1
    u = fourier_approx(d[1], H, L; cc=cc, N=N)
    F = SteadyWaves.wave_power(u, N) / √k^5 # F / ρ√g³
    T = wave_period(u, d[1], N) * √g # T * √g
    for i in eachindex(d)
        if i>1
            fourier_approx!(u, d[i], d[i-1], F / √d[i]^5, T / √d[i];  cc=cc, N=N)
            K[i] = u[2N+7] / u[2N+3] * d[i] / H
        end
    end
    return K
end

"""
    fourier_approx!(u, d, d_p, F, T; cc=1, N=10, g=9.81)

Update approximate solution `u` of a steady wave of power `F` and period `T`
propagating in water of changing depth from `d` to `d_p` using Fourier Approximation Method.

# Arguments
- `u`: solution matrix being mutated
- `d`: initial water depth (m)
- `d_p`: target water depth (m)
- `F`: wave power (kg m/s)
- `T`: wave period (s)
- `cc`: current criterion; `cc=1`, `cc=CC_STOKES` - Stokes (default), `cc=2`, `cc=CC_EULER` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
"""
function fourier_approx!(u, d, d_p, F, T; cc=CC_STOKES, N=10)
    init_conditions!(d_p / d, u, N)
    params = [F, T]

    nonlinear_system(du,u,p) = nonlinear_system_shoaling(
        du, u, p,
        current_criterion(cc)
    )

    problem = NonlinearProblem(nonlinear_system, u[1:2N+7], params)
    solution = solve(problem, RobustMultiNewton())
    u[1:2N+7] = solution.u
    return nothing
end

function init_conditions!(ratio_d, u, N)
    u[1:N+1] =  1 .+ (u[1:N+1] .- 1) / ratio_d # kη
    u[N+2:2N+1] /= √ratio_d # B
    u[2N+2] /= √ratio_d # c√(k/g)
    u[2N+3] /= ratio_d # kη̄
    u[2N+4] /= √ratio_d^3 # q√(k³/g)
    u[2N+5] /= ratio_d # rk/g
    u[2N+6] /= √ratio_d # Ū√(k/g)
    u[2N+7] /= ratio_d # kH
    return nothing
end

function wave_height_condition(u,N)
    return u[1] - u[N+1] - u[2N+7]
end
function power_condition(u,N,p)
    return wave_power(u, N) - p * √u[2N+3]^5
end

"""
    nonlinear_system_shoaling(du, u, p)

Define nonlinear system for shoaling waves `f(u) = 0` with parameters `p`.

"""
function nonlinear_system_shoaling(du, u, p, cc_equation)
    N = (length(u) - 7) ÷ 2
    
    for m in 0:N
        du[m+1] = kinematic_surface_condition(u, N, m)
        du[N+1+m+1] = dynamic_surface_condition(u, N, m)
    end

    du[2N+3] = mean_depth(u,N)
    
    du[2N+4] = wave_height_condition(u,N)

    du[2N+5] = period_condition(u,N,p[2])

    du[2N+6] = cc_equation(u,N)

    du[2N+7] = power_condition(u,N,p[1])
    return nothing
end
