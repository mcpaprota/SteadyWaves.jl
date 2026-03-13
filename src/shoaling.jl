# SPDX-License-Identifier: MIT

# Functions for shoaling calculations based on Fourier Approximation Method
module Shoaling

using ..Output
using ..Output: wave_power, wave_period
using ..Params
using ..Physics
using NonlinearSolve
using ..Steady: fourier_approx
using ..NonlinearSystem: nonlinear_system_base!, period_condition, power_condition, current_condition_factory, height_condition,
        kinematic_surface_condition, dynamic_surface_condition, mean_depth_condition,
        ConditionStruct
"""
    topo_approx(d, H, L; cc=2, N=10, g=G)

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
function topo_approx(d, H, L; cc=CC_STOKES, N=10, g=G)
    k = 2π / L # initial wave number (rad/m
    ω = √(g * k * tanh(k * d[1])) # initial angular wave frequency (rad/s)

    K = zero(float(d))
    K[1] = 1
    u = fourier_approx(d[1], H, L; cc=cc, N=N)
    F = wave_power(u, N) / √k^5 # F / ρ√g³
    T = wave_period(u, d[1], N) * √g # T * √g
    for i in eachindex(d)
        if i>1
            fourier_approx!(u, d[i], d[i-1], F / √d[i]^5, T / √d[i];  cc=cc, N=N)
            K[i] = u[2N+H_INDEX] / u[2N+D_INDEX] * d[i] / H
        end
    end
    return K
end

"""
    fourier_approx!(u, d, d_p, F, T; cc=1, N=10, g=G)

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


    eta = Output.Elevation.DirectElevation

    _kinematic_surface_condition(u,N,m) = kinematic_surface_condition(u,N,m,eta)
    _dynamic_surface_condition(u,N,m) = dynamic_surface_condition(u,N,m,eta)
    _mean_depth_condition(u,N) = mean_depth_condition(u,N,eta)

    _height_condition(u,N) = height_condition(u,N,eta)

    _period_condition(u,N) = period_condition(u, N, T)
    _power_condition(u,N) = power_condition(u,N,F)
    _current_condition = current_condition_factory(cc)

    conditions = [
        ConditionStruct(_kinematic_surface_condition, 0:N),
        ConditionStruct(_dynamic_surface_condition, 0:N),
        ConditionStruct(_mean_depth_condition),
        ConditionStruct(_period_condition),
        ConditionStruct(_current_condition),
        ConditionStruct(_height_condition),
        ConditionStruct(_power_condition)
    ]

    _nonlinear_system!(du,u,p) = nonlinear_system_base!( du, u, N, conditions)

    problem = NonlinearProblem(_nonlinear_system!, u[1:2N+7])
    solution = solve(problem, RobustMultiNewton())
    u[1:2N+7] = solution.u
    return nothing
end

function init_conditions!(ratio_d, u, N)
    u[elevation_indexes(N)] =  1 .+ (u[elevation_indexes(N)] .- 1) / ratio_d # kη
    u[stream_indexes(N)] /= √ratio_d # B
    u[2N+C_INDEX] /= √ratio_d # c√(k/g)
    u[2N+D_INDEX] /= ratio_d # kη̄
    u[2N+Q_INDEX] /= √ratio_d^3 # q√(k³/g)
    u[2N+R_INDEX] /= ratio_d # rk/g
    u[2N+U_INDEX] /= √ratio_d # Ū√(k/g)
    u[2N+H_INDEX] /= ratio_d # kH
    return nothing
end

end