# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)
module Steady

using ..Index
using ..Surface
using ..Wave: WaveStruct, Wave
using ..Output
using ..Params
using ..Physics
using ..NonlinearSystem: fourier_approx_base, parameter_condition_constant, parameter_condition_factory,
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
function fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=G,rho=RHO)
    L , T = Params.L(P,pc), Params.T(P,pc)

    idx = Index.default_indexes(N)

    # create default compiler
    compiler = WaveStruct(idx)

    # create dimensional factor compiler 
    df_compiler = Wave.dimensional_factor_compiler(compiler.D, d, g, rho)

    # set dimensionless height, length and period from dimensional value
    # if property is nothing given change is skipped
    compiler = WaveStruct(compiler, df_compiler;
        H = H/M,
        L = L,
        T = T,
    )

    u = init_conditions(d, H, P, Int(pc), N, M, idx)

    conditions = [
        ConditionStruct(kinematic_surface_condition,0:N),
        ConditionStruct(dynamic_surface_condition,0:N),
        ConditionStruct(mean_depth_condition),
        ConditionStruct(parameter_condition_factory(pc)),
        ConditionStruct(current_condition_factory(cc)),
        ConditionStruct(height_condition)
    ]


    for m in 1:M
        # update height
        compiler = WaveStruct(compiler, df_compiler; H = H * m/M)

        w = fourier_approx_base(u,compiler,conditions)
        u[:] = w.raw
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
function init_conditions(d, H, P, pc, N, M, idx)
    return dimensionless_linear_solution(d,H / M, P, pc, idx)
end

function init(d,P,pc,idx)
    k = Int(pc) == Int(PC_LENGTH) ? 2π / P : linear_wave_number(d, 2π / P) # wave number (rad/s)
    u0 = zeros(idx.U)

    return k, u0
end

function dimensionless_linear_solution(d, H, P, pc, idx)
    k, u0 = init(d, P, pc, idx)

    u0[idx.eta] = @. k * d + 0.5 * k * H * cos((0:idx.N) * π / idx.N) # kη
    u0[idx.psi[begin]] = 0.5 * k * H / √tanh(k * d) # Bk/g
    u0[idx.C] = √tanh(k * d) # c√(k/g)
    u0[idx.D] = k * d # kη̄
    u0[idx.Q] = 0 # q√(k³/g)
    u0[idx.R] = tanh(k * d) / 2 # rk/g
    u0[idx.U] = √tanh(k * d) # Ū√(k/g)
    return u0
end

function dimensional_linear_solution(d, H, P, pc, idx, g = G)
    k, u0 = init(d, P, pc, idx)

    u0[idx.eta] = @.d + 0.5 * H * cos((0:N) * π / idx.N) # η
    u0[idx.psi[begin]] = 0.5 * g * H / √tanh(k * d) # B 

    u0[idx.C] = √(g *tanh(k * d) / k) # c
    u0[idx.D] = d # η̄
    u0[idx.Q] = 0 # q
    u0[idx.R] = 0.5 * g * tanh(k * d) / k # r 
    u0[idx.U] = √(g * tanh(k * d) / k) # Ū
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