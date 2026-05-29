# SPDX-License-Identifier: MIT

# Basic functions for Fourier Approximation Method (FAM)
module Steady

using ..Index
using ..Surface
using ..Wave: WaveStruct, Wave
using ..DimensionalFactor: dimensional_factor_compiler
using ..Output
using ..Indirect
using ..Params
using ..Physics
using ..Linear
using ..Current
using ..NonlinearSystem: fourier_approx_base, ConditionStruct
using ..Condition: parameter_condition_factory,
    current_condition_factory, height_condition,
    kinematic_surface_condition,
    mean_depth_condition, dynamic_condition_factory,condition_factory
"""
    fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=G,rho=RHO,sigma=SIGMA,
        eta_type::ElevationType = Params.FOURIER_ELEVATION,
        wave_type::Params.WaveType = Params.GRAVITY_WAVE,
        deep_water = false,
        c_e::Union{Nothing,Number}=nothing,
    )
Approximate solution `w` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using Fourier Approximation Method.

# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=PC_LENGTH` - length (default), `pc=PC_PERIOD` - period
- `cc`: current criterion; `cc=CC_STOKES` - Stokes (default), `cc=CC_EULER` - Euler
- `N`: number of solution eigenvalues, defaults to `N=10`
- `M`: number of height steps, defaults to `M=1`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`
- `rho`: density (kg/m^3),default to `rho=1000`
- `sigma` - surface tension coefficient `sigma=0.073`
- `eta_type` - elevation type - FOURIER_ELEVATION, DIRECT_ELEVATION
- `deep_water`: depth flag - `true` - infinite depth `false` finite depth
- `c_e`: eulerian current - if `cc = CC_ARBITRARY` describe eulerian current
# Output
- `w`: Wave Structure
- `df`: Dimensional Factor Structure
"""
function fourier_approx(d, H, P; pc=PC_LENGTH, cc=CC_STOKES, N=10, M=1, g=G,rho=RHO,sigma=SIGMA,
    eta_type::ElevationType = Params.FOURIER_ELEVATION,
    wave_type::Params.WaveType = Params.GRAVITY_WAVE,
    deep_water = false,
    c_e::Union{Nothing,Number}=nothing,
    )

    config = Params.ConfigStruct(
        cc=cc,
        eta_type=eta_type,
        deep_water=deep_water,
        wave_type=wave_type
    )

    c_e = config.cc == Params.CC_ARBITRARY ? nothing : c_e

    L , T = Params.L(P,pc), Params.T(P,pc)

    definition = Params.Definition(
        d = d,
        H = H,
        L = L,
        T = T,
        c_e = c_e
    )

    

    physics = Physics.PhysicsStruct(g,rho,sigma) 


    validate_config(d,H,L,T,config,physics,N,M)

    return fourier_approx(definition,config, physics; N=N,M=M)
end

function validate_config(d,H,L,T,config,physics,N,M)
    Physics.validate_constants(physics)

    Physics.validate_parameters(H,L,T,d)

    @assert (L === nothing) || (T === nothing)

    @assert N > 1

    @assert M > 0
end

function fourier_approx(definition::Params.Definition, config::Params.ConfigStruct, physics::Physics.PhysicsStruct; N=10, M=1)
    idx = Index.default_indexes(N)

    # Calculate initial wave using linear theory
    def_0 = Wave.set_values(definition,Params.Definition(H=definition.H/M))
    w, df = Linear.linear_solution(def_0, config, physics, idx)

    # Approximate nonlinear solution then increase wave height up to selected
    for m in 1:M
        # update height
        def_m = Wave.set_values(definition,Params.Definition(H=definition.H*m/M))

        w, df = conformal_furrier_approx(def_m,config,physics,idx,w)
    end


    w = output_wave(w,idx,config)

    push!(w.raw,w.H)
    return w, df

end

function dimensionless_fourier_approx(definition::Params.Definition,config::Params.ConfigStruct;dimensionless_sigma=0, N=10)
    idx = Index.dynamic_indexes(N,
        eta=1:N+1,
        psi=1:N,
        C= config.indirect_celerity == false,
        R=true,
        U=true,
        Q=true
    )

    # create default compiler
    compiler = WaveStruct(idx,config)

    # set dimensionless height and sigma
    compiler = Wave.set_compilator_values(
        compiler,
        WaveStruct(
            D = definition.d,
            H = definition.H,
            sigma = dimensionless_sigma,
            L = 2pi,
        )
    )

    # initial conditions
    w, _ = Linear.dimensionless_linear_solution(config, idx, compiler)

    conditions = [
        ConditionStruct(kinematic_surface_condition,0:N),
        ConditionStruct(dynamic_condition_factory(config),0:N),
        ConditionStruct(current_condition_factory(config)),
        ConditionStruct(mean_depth_condition),
        ConditionStruct(height_condition)
    ]

    conditions  = filter(x -> x.condition !== nothing, conditions)

    println.(conditions)

    w = fourier_approx_base(w.raw,compiler,conditions)

    w = output_wave(w,idx,config)

    return w, nothing

end

function output_wave(w,idx,config)
    if config.indirect_celerity
        w = Wave.set_values(w,WaveStruct(C = Indirect.indirect_celerity_factory(config)(w)))
    end

    return Wave.set_values(w,
        WaveStruct(
            eta = Surface.struct_with_derived_values(w.eta,idx,config.eta_type),
            P = (kx,kz) -> Indirect.indirect_pressure(w,kx,kz),
            F = Indirect.indirect_wave_power(w),
            L = w.L === nothing ? Indirect.indirect_wavelength(w) : nothing,
            T = w.T === nothing ? Indirect.indirect_wave_period(w) : nothing,
        )
    )
end
# update wave to fit definition
function conformal_furrier_approx(definition,config,physics,idx,w)

    compiler = WaveStruct(idx,config)

    # create dimensional factor compiler 
    df_compiler = dimensional_factor_compiler(definition.d, physics)

    # set dimensionless height, length and period from dimensional value
    compiler = Wave.set_compilator_values(
        compiler,
        WaveStruct(definition,physics),
        df_compiler,
    )

    conditions = condition_factory(definition,config,idx.N)

    w = fourier_approx_base(w.raw,compiler,conditions)

    w = output_wave(w,idx,config)

    return w, WaveStruct(w.raw, df_compiler, compiler)

end

end