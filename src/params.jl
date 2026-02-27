module Params

@enum CurrentCriterion begin
    CC_STOKES = 1
    CC_EULER = 2
end

@enum ParameterCriterion begin
    PC_LENGTH = 1
    PC_PERIOD= 2
end
export CurrentCriterion, CC_EULER, CC_STOKES

export ParameterCriterion, PC_LENGTH, PC_PERIOD
end