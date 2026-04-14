module Params

@enum CurrentCriterion begin
    CC_STOKES = 1
    CC_EULER = 2
    CC_INVALID = 0
end

@enum ParameterCriterion begin
    PC_LENGTH = 1
    PC_PERIOD = 2
    PC_INVALID = 0
end

@enum ElevationType begin
    DIRECT_ELEVATION = 1
    FOURIER_ELEVATION = 2
    INVALID_ELEVATION = 0
end

@enum WaveType begin
    GRAVITY_WAVE = 1
    CAPILLARY_WAVE = 2
    GRAVITY_CAPILLARY_WAVE = 3
    INVALID_WAVE = 0
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
    deep_water::Bool

    ConfigStruct(;
        cc=CC_INVALID,
        pc=PC_INVALID,
        eta_type=INVALID_ELEVATION,
        wave_type=GRAVITY_WAVE,
        deep_water::Bool=false
    ) = new(
        typeof(cc) == Int ? CurrentCriterion(cc) : cc,
        typeof(pc) == Int ? ParameterCriterion(pc) : pc,
        typeof(eta_type) == Int ? ElevationType(eta_type) : eta_type,
        typeof(wave_type) == Int ? ElevationType(wave_type) : wave_type,
        deep_water,
    )
end

export CurrentCriterion, CC_EULER, CC_STOKES

export ParameterCriterion, PC_LENGTH, PC_PERIOD

export ElevationType

export WaveType

end