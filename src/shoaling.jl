# SPDX-License-Identifier: MIT

# Functions for shoaling calculations based on Fourier Approximation Method
module Shoaling

using ..Index
using ..Surface
using ..Wave
using ..Output
using ..Output: wave_power, wave_period
using ..Params
using ..DimensionalFactor: dimensional, dimensional_factor_compiler
using ..Physics
using ..Wave: WaveStruct
using ..Steady: fourier_approx
using ..NonlinearSystem: fourier_approx_base, ConditionStruct
using ..Condition: period_condition, power_condition, current_condition_factory, height_condition,
        kinematic_surface_condition, mean_depth_condition,
        dynamic_condition_factory
"""
    topo_approx(d, H, L; cc=CC_STOKES, N=10, g=G, rho = RHO, sigma=0,eta_type=Params.FOURIER_ELEVATION,c_e=nothing)

Calculate shoaling coefficients `K` in range of depth values `d`
for wave of length `L` and height `H`.

# Arguments
- `d`: vector of decreasing water depths (m)
- `L`: initial wavelength (m) - corresponding to d[1]
- `H`: initial wave height (m) - corresponding to d[1]
- `cc`: current criterion; `cc=CC_STOKES` - Stokes (default), `cc=CC_EULER` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`
- `rho`: density (kg/m^3),default to `rho=1000`
- `sigma` - surface tension coefficient `sigma=0.073`
- `eta_type` - elevation type - FOURIER_ELEVATION, DIRECT_ELEVATION
- `deep_water`: depth flag - `true` - infinite depth `false` finite depth
- `c_e`: eulerian current - if `cc = CC_ARBITRARY` describe eulerian current

# Output
- `K`: vector of shoaling coefficient values
"""
function topo_approx(d, H, L; cc=CC_STOKES, N=10, g=G, rho = RHO, sigma=0,
        eta_type=Params.FOURIER_ELEVATION,
        c_e=nothing,
    )
    
    config = Params.ConfigStruct(pc=PC_LENGTH, cc=cc, eta_type=eta_type)

    definition = Params.Definition(
        d = d[1],
        H = H,
        L = L,
        c_e = c_e
    )
   
    physics = Physics.PhysicsStruct(g,rho,sigma)

    return topo_approx(definition,config,physics,d,N=N)
end

function topo_approx(definition, config,physics,d; N=10)
    H = definition.H
    idx = Index.default_indexes(N)
    k = 2π / definition.L # initial wave number (rad/m)

    K = zero(float(d))
    K[1] = 1
    w, df = fourier_approx(definition,config,physics; N=N)

    wd = dimensional(w,df)

    for i in eachindex(d)[begin+1:end]
        w, df = update_depth_fourier_approx(w, d[i], d[i-1], wd.F, wd.T, idx, config, physics; N=N)
        K[i] = w.H / df.H / H
    end
    return K
end
"""
    update_depth_fourier_approx(w, d, d_p, F, T, idx; cc=CC_STOKES, N=10, g=G,rho=RHO,eta_type = Params.FOURIER_ELEVATION)

Approximate solution `w` of a steady wave that approaches change in depth

# Arguments
- `w`: previous wave
- `d`: water depth (m)
- `d_p`: previous water depth (m)
- `F`: wave flux
- `T`: wave period (s)
- `idx`: solution array indexes
- `cc`: current criterion; `cc=CC_STOKES` - Stokes (default), `cc=CC_EULER` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`
- `rho`: density (kg/m^3),default to `rho=1000`
- `sigma` - surface tension coefficient `sigma=0.073`
- `eta_type` - elevation type - FOURIER_ELEVATION, DIRECT_ELEVATION

# Output
- `w`: Wave Structure
- `df`: Dimensional Factor Structure
"""
function update_depth_fourier_approx(w, d, d_p, F, T, idx; cc=CC_STOKES, N=10, g=G,rho=RHO,eta_type = Params.FOURIER_ELEVATION)
    
    config = Params.ConditionStruct(cc=cc, eta_type=eta_type)

    physics = Physics.PhysicsStruct(g,rho,sigma)

    return update_depth_fourier_approx(w,d,d_p,F,T,idx,config,physics,N=N)
end

function update_depth_fourier_approx(w, d, d_p, F, T, idx,config,physics; N=10)
    init_conditions!(d_p / d, w.raw, idx,config)

    # create default compiler
    compiler = Wave.WaveStruct(idx,config)

    # create dimensional_factor_compiler 
    df_compiler = dimensional_factor_compiler(d, physics)

    # set dimensionless period and wave_power from dimensional values with respect to depth
    compiler = Wave.set_compilator_values(compiler,
        Wave.WaveStruct(
            T = T,
            F = F,
        ),
        df_compiler
    )

    conditions = [
        ConditionStruct(kinematic_surface_condition, 0:N),
        ConditionStruct(dynamic_condition_factory(config), 0:N),
        ConditionStruct(mean_depth_condition),
        ConditionStruct(period_condition),
        ConditionStruct(current_condition_factory(config)),
        ConditionStruct(height_condition),
        ConditionStruct(power_condition)
    ]

    w = fourier_approx_base(w.raw,compiler,conditions)

    df = WaveStruct(w.raw,df_compiler, compiler)

    return w, df
end

function init_conditions!(ratio_d, u, idx,config)

    if config.eta_type == Params.DIRECT_ELEVATION
        
        u[idx.eta] =  1 .+ (u[idx.eta] .- 1) / ratio_d # kη
        
    elseif config.eta_type == Params.FOURIER_ELEVATION

        u[idx.eta[begin]] =  1 .+ (u[idx.eta[begin]] .- 1) / ratio_d # kη
        u[idx.eta[begin+1:end]] /= ratio_d

    end
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