# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)
module Steady

using ..Index
using ..Surface
using ..Wave: WaveStruct, Wave
using ..DimensionalFactor: dimensional_factor_compiler
using ..Output
using ..Params
using ..Physics
using ..Linear
using ..NonlinearSystem: fourier_approx_base, ConditionStruct
using ..Condition: parameter_condition_factory,
    current_condition_factory, height_condition,
    kinematic_surface_condition,
    mean_depth_condition, dynamic_condition_factory
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
function fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=G,rho=RHO,sigma=SIGMA,
    eta_type::ElevationType = Params.FOURIER_ELEVATION,
    deep_water::Bool = false
    )

    config = Params.ConfigStruct(
        cc=cc, pc=pc,
        eta_type=eta_type,
        deep_water=deep_water,
    )

    physics = Physics.PhysicsStruct(g,rho,sigma) 

    return fourier_approx(d,H,P,config, physics; N=N,M=M)
end


function fourier_approx(d, H, P,config::Params.ConfigStruct, physics::Physics.PhysicsStruct; N=10, M=1)

    L , T = Params.L(P,config.pc), Params.T(P,config.pc)

    idx = Index.default_indexes(N)

    # create default compiler
    compiler = WaveStruct(idx,config)

    # create dimensional factor compiler 
    df_compiler = dimensional_factor_compiler(d, physics)

    # set dimensionless height, length and period from dimensional value
    # if property is nothing given change is skipped
    compiler = WaveStruct(compiler, df_compiler;
        H = H/M,
        L = L,
        T = T,
        sigma = physics.sigma
    )

    # initial conditions
    w, _ = Linear.linear_solution(d, P, config, idx, compiler, df_compiler)

    conditions = [
        ConditionStruct(kinematic_surface_condition,0:N),
        ConditionStruct(dynamic_condition_factory(config),0:N),
        ConditionStruct(mean_depth_condition),
        ConditionStruct(parameter_condition_factory(config.pc)),
        ConditionStruct(current_condition_factory(config.cc)),
        ConditionStruct(height_condition)
    ]

    for m in 1:M
        # update height
        compiler = WaveStruct(compiler, df_compiler; H = H * m/M)

        w = fourier_approx_base(w.raw,compiler,conditions)
    end


    w = WaveStruct(w;
        eta = Surface.struct_with_derived_values(w.eta,idx,config.eta_type)
    )

    push!(w.raw,w.H)
    return w, WaveStruct(w.raw, df_compiler, compiler)

end

end