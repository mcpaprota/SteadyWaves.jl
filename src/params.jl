module Params

@enum CurrentCriterion begin
    CC_STOKES = 1
    CC_EULER = 2
end

@enum ParameterCriterion begin
    PC_LENGTH = 1
    PC_PERIOD = 2
end

@enum ElevationType begin
    DIRECT_ELEVATION = 1
    FOURIER_ELEVATION = 2
end

function T(P,pc)
    return Int(pc) == Int(PC_PERIOD) ? P : nothing
end

function L(P,pc)
    return Int(pc) == Int(PC_LENGTH) ? P : nothing
end

export CurrentCriterion, CC_EULER, CC_STOKES

export ParameterCriterion, PC_LENGTH, PC_PERIOD

export ElevationType
end