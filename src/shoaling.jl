# SPDX-License-Identifier: MIT

# Functions for shoaling calculations based on Fourier Approximation Method
module Shoaling

using ..Index
using ..Surface
using ..Wave
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
    idx = Index.default_indexes(N)
    k = 2π / L # initial wave number (rad/m
    ω = √(g * k * tanh(k * d[1])) # initial angular wave frequency (rad/s)

    K = zero(float(d))
    K[1] = 1
    u = fourier_approx(d[1], H, L; cc=cc, N=N)
    F = wave_power(u, N) / √k^5 # F / ρ√g³
    T = wave_period(u, d[1], N) * √g # T * √g
    for i in eachindex(d)
        if i>1
            fourier_approx!(u, d[i], d[i-1], F, T, idx;  cc=cc, N=N)
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
function fourier_approx!(u, d, d_p, F, T, idx; cc=CC_STOKES, N=10)
    init_conditions!(d_p / d, u, idx)

    # create default compiler
    compiler = Wave.WaveStruct(idx)

    # create dimensional_factor_compiler 
    df_compiler = Wave.dimensional_factor_compiler( d, 1, 1)

    # set dimensionless period and wave_power from dimensional values with respect to depth
    compiler = Wave.WaveStruct(compiler, df_compiler;
        T = T,
        F = F,
    )

    conditions = [
        ConditionStruct(kinematic_surface_condition, 0:N),
        ConditionStruct(dynamic_surface_condition, 0:N),
        ConditionStruct(mean_depth_condition),
        ConditionStruct(period_condition),
        ConditionStruct(current_condition_factory(cc)),
        ConditionStruct(height_condition),
        ConditionStruct(power_condition)
    ]

    _nonlinear_system!(du,u,p) = nonlinear_system_base!( du, u, conditions, compiler)

    problem = NonlinearProblem(_nonlinear_system!, u[1:2N+7])
    solution = solve(problem, RobustMultiNewton())
    u[1:2N+7] = solution.u
    return nothing
end

function init_conditions!(ratio_d, u, idx)
    u[idx.eta] =  1 .+ (u[idx.eta] .- 1) / ratio_d # kη
    u[idx.psi] /= √ratio_d # B
    u[idx.C] /= √ratio_d # c√(k/g)
    u[idx.D] /= ratio_d # kη̄
    u[idx.Q] /= √ratio_d^3 # q√(k³/g)
    u[idx.R] /= ratio_d # rk/g
    u[idx.U] /= √ratio_d # Ū√(k/g)
    u[idx.H] /= ratio_d # kH
    return nothing
end

end