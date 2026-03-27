module Condition

using ..Wave:WaveStruct
using ..Output
using ..Output: indirect_wave_power, surface_stream_eigenfunction, dimensionless_pressure
using ..Params


function mean_depth_condition(w::WaveStruct)
    return w.eta.avg - w.D
end

function kinematic_surface_condition(w::WaveStruct,m)
    psi = sum([surface_stream_eigenfunction(sinh, cos, w.raw, w.N, m, j) for j in 1:w.N])
    return psi - w.U * (w.eta.point(m) - w.D) - w.Q
end

function dynamic_surface_condition(w::WaveStruct,m)
    kx = m/w.N * pi
    kz = w.eta.point(m)

    return dimensionless_pressure(w.raw, w.N, kx, kz)
end

function height_condition(w::WaveStruct, p)
    return w.eta.max - w.eta.min - w.D * p
end

function height_condition(w::WaveStruct)
    return w.eta.max - w.eta.min - w.H
end

function power_condition(w::WaveStruct)
    return indirect_wave_power(w) - w.F
end

function euler_condition(w::WaveStruct)
    return w.U - w.C
end

function stokes_condition(w::WaveStruct)
    return euler_condition(w::WaveStruct) - w.Q / w.D
end

function length_condition(w::WaveStruct)
    return w.L - 2π 
end

function period_condition(w::WaveStruct)
    return w.C * w.T- 2π   
end


function current_condition_factory(cc)
    if Int(cc) == Int(CC_STOKES)
        return stokes_condition
    elseif Int(cc) == Int(CC_EULER)
        return euler_condition
    else
        throw(error("Unknown current criterion $cc"))
    end
end

function parameter_condition_factory(pc)
    if Int(pc) == Int(PC_LENGTH)
        return length_condition
    elseif Int(pc) == Int(PC_PERIOD)
        return period_condition
    else
        throw(error("Unknown parameter criterion $pc"))
    end
end

end