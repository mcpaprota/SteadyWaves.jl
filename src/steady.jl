# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)
module Steady

using ..Output
using ..Params
using ..Physics
using NonlinearSolve
using ..NonlinearSystem: nonlinear_system_base!, parameter_condition_constant, parameter_condition_factory,
    current_condition_factory, height_condition,
    kinematic_surface_condition, dynamic_surface_condition, mean_depth_condition,
    ConditionStruct
"""
    fourier_approx(d, H, P; pc=1, cc=1, N=10, M=1, g=G)

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
function fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=G)
    u = init_conditions(d, H, P, Int(pc), N, M)

    parameter_constant = parameter_condition_constant(pc, P, d, g)
    _parameter_condition(u,N) = parameter_condition_factory(pc)(
        u,N,
        parameter_constant
    )
    _current_condition = current_condition_factory(cc)

    for m in 1:M
        _height_condition(u, N) = height_condition(u, N, H / d * m / M)

        conditions = [
            ConditionStruct(kinematic_surface_condition,0:N),
            ConditionStruct(dynamic_surface_condition,0:N),
            ConditionStruct(mean_depth_condition),
            ConditionStruct(_parameter_condition),
            ConditionStruct(_current_condition),
            ConditionStruct(_height_condition)
        ]

        _nonlinear_system!(du,u,p) = nonlinear_system_base!(du,u,N,conditions)

        problem = NonlinearProblem(_nonlinear_system!, u)
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
    Int(pc) == Int(PC_LENGTH) ? k = 2π / P : k = linear_wave_number(d, 2π / P) # wave number (rad/s)
    u0 = zeros(2N + U_INDEX)
    u0[elevation_indexes(N)] = @. k * d + 1 / 2 * k * H / M * cos((0:N) * π / N) # kη
    u0[stream_indexes(N)] = [k * H / M / 2 / √tanh(k * d); zeros(N - 1)] # B
    u0[2N+C_INDEX] = √tanh(k * d) # c√(k/g)
    u0[2N+D_INDEX] = k * d # kη̄
    u0[2N+Q_INDEX] = 0 # q√(k³/g)
    u0[2N+R_INDEX] = tanh(k * d) / 2 # rk/g
    u0[2N+U_INDEX] = √tanh(k * d) # Ū√(k/g)
    return u0
end

"""

    linear_wave_number(d, ω, g=G, ϵ=10^-12)

Calculate linear_wave_number `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""

function linear_wave_number(d, ω, g=G, ϵ=10^-12)
    k = k₀ = ω^2 / g # initial guess
    while max(abs(k * tanh(k * d) - k₀)) > ϵ
        k = k₀ / tanh(k * d)
    end
    return k
end

end