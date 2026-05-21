module Params

"""
Current Criterion selects condition that relates `c` to `u`

- CC_EULER - euler condition `u - c = 0`
- CC_STOKES - stokes condition `u - c - Q/d = 0`
- CC_ARBITRARY - current condition `u -c - c_e = 0`
- CC_INVALID - for internal use only. Most likely value was:
    - not selected
    - used in invalid context
"""
@enum CurrentCriterion begin
    CC_STOKES = 1
    CC_EULER = 2
    CC_ARBITRARY = 3
    CC_INVALID = 0
end
"""
Parameter Criterion selects condition that enforce wave length via `L` or `T`

- PC_LENGTH - length condition `L - 2π = 0`
- PC_PERIOD - period condition `T*c - 2π = 0`
- PC_INVALID - for internal use only. Most likely value was:
    - not selected
    - used in invalid context
"""
@enum ParameterCriterion begin
    PC_LENGTH = 1
    PC_PERIOD = 2
    PC_INVALID = 0
end
"""
Elevation Type selects representesion of free surface `η`

- DIRECT_ELEVATION - represents `η` as `N+1` points along a wave
- FOURIER_ELEVATION - represents `η` as `N+1` cosine frequencies of a wave
- INVALID_ELEVATION - for internal use only. Most likely value was:
    - not selected
    - used in invalid context
"""
@enum ElevationType begin
    DIRECT_ELEVATION = 1
    FOURIER_ELEVATION = 2
    INVALID_ELEVATION = 0
end
"""
Wave Type selects dominant forces shaping wave

- GRAVITY_WAVE              - apply only gravitational force to a wave

- CAPILLARY_WAVE            - apply only capillary force to a wave 

- GRAVITY_CAPILLARY_WAVE    - apply both capilary and gravitational forces to a wave

- INVALID_WAVE - for internal use only. Most likely value was:
    - not selected
    - used in invalid context
"""
@enum WaveType begin
    GRAVITY_WAVE = 1
    CAPILLARY_WAVE = 2
    GRAVITY_CAPILLARY_WAVE = 3
    INVALID_WAVE = 0
end
"""
DepthType selects stream function representesion
- DEEP_WATER - selects representesion stable in deep water - `e^y`
    
- SHALLOW_WATER - selected representesion stable in shallow water `sinh(y+d)/cosh(d)`
- STABLE - selects representesion stable in all conditions but slow `sinh(y)+tanh(d)*cosh(y)`
"""
@enum DepthType begin
    DEEP_WATER = 1
    SHALLOW_WATER = 2
    STABLE = 3
end

function parse_depth_type(deep::Bool)
    return deep ? DEEP_WATER : SHALLOW_WATER
end

function T(P,pc)
    return Int(pc) == Int(PC_PERIOD) ? P : nothing
end

function L(P,pc)
    return Int(pc) == Int(PC_LENGTH) ? P : nothing
end

struct ConfigStruct
    cc::CurrentCriterion
    pc::ParameterCriterion
    eta_type::ElevationType
    wave_type::WaveType
    deep_water::DepthType
    indirect_celerity::Bool

    ConfigStruct(;
        cc=CC_INVALID,
        pc=PC_INVALID,
        eta_type=INVALID_ELEVATION,
        wave_type=GRAVITY_WAVE,
        deep_water=SHALLOW_WATER,
        indirect_celerity::Bool=false,
    ) = new(
        typeof(cc) == Int ? CurrentCriterion(cc) : cc,
        typeof(pc) == Int ? ParameterCriterion(pc) : pc,
        typeof(eta_type) == Int ? ElevationType(eta_type) : eta_type,
        typeof(wave_type) == Int ? ElevationType(wave_type) : wave_type,
        typeof(deep_water) == Bool ? parse_depth_type(deep_water) : deep_water,
        indirect_celerity,
    )
end

struct Definition
    d
    H
    F
    L
    T
    c_e

    Definition(;
        d=nothing,
        H=nothing,
        F=nothing,
        L=nothing,
        T=nothing,
        c_e=nothing
    ) =new(d, H, F, L, T, c_e)
end

export CurrentCriterion, CC_EULER, CC_STOKES

export ParameterCriterion, PC_LENGTH, PC_PERIOD

export ElevationType

export WaveType

end