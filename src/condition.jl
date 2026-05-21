module Condition

using ..Wave:WaveStruct
using ..Output
using ..Indirect: indirect_wave_power, indirect_pressure, indirect_surface_tension,indirect_dynamic_pressure
using ..Params


function mean_depth_condition(w::WaveStruct)
    return w.eta.avg - w.D
end

function still_water_depth_condition(w::WaveStruct)
    return w.Q - w.U * (w.eta.avg - w.D)
end

function depth_condition_factory(config)
    dict = Dict{Params.ReferenceLevel,Function}(
        Params.MEAN_DEPTH => mean_depth_condition,
        Params.STILL_DEPTH => still_water_depth_condition,
    )
    return dict[config.reference_level]
end

function kinematic_surface_condition(w::WaveStruct,m)
    kx = w.eta.point.x(m)
    kz = w.eta.point(m)

    return w.v.psi(kx,kz) - w.U * (kz - w.D) - w.Q
end


function gravity_dynamic_surface_condition(w::WaveStruct,m)
    kx = w.eta.point.x(m)
    kz = w.eta.point(m)

    return indirect_pressure(w,kx,kz)
end

function gravity_capillary_dynamic_surface_condition(w::WaveStruct,m)
    kx = w.eta.point.x(m)
    kz = w.eta.point(m)

    dz_dx_1 = w.eta.point.dz_dx_1(m)
    dz_dx_2 = w.eta.point.dz_dx_2(m)

    return indirect_pressure(w,kx,kz) - indirect_surface_tension(w.sigma, dz_dx_1, dz_dx_2)
end

function capillary_dynamic_surface_condition(w::WaveStruct,m)
    kx = m/w.N * pi
    kz = w.eta.point(m)

    dz_dx_1 = w.eta.point.dz_dx_1(m)
    dz_dx_2 = w.eta.point.dz_dx_2(m)

    return indirect_dynamic_pressure(w,kx,kz) - indirect_surface_tension(w.sigma, dz_dx_1, dz_dx_2)
end


function dynamic_condition_factory(config::Params.ConfigStruct)
    if config.wave_type == Params.GRAVITY_WAVE

        return gravity_dynamic_surface_condition

    elseif config.wave_type == Params.GRAVITY_CAPILLARY_WAVE

        return gravity_capillary_dynamic_surface_condition

    elseif config.wave_type == Params.CAPILLARY_WAVE

        return capillary_dynamic_surface_condition

    else
        throw(error("Unknown wave type $(config.wave_type)"))
    end

end

function height_condition(w::WaveStruct)
    return w.eta.max - w.eta.min - w.H
end

function power_condition(w::WaveStruct)
    return indirect_wave_power(w) - w.F
end

function length_condition(w::WaveStruct)
    return w.L - 2π 
end

function period_condition(w::WaveStruct)
    return w.C * w.T- 2π   
end

function current_condition(w::WaveStruct)
    return w.C - w.U - w.c_e
end

function current_condition_factory(config::Params.ConfigStruct)
    if config.indirect_celerity
        return nothing
    else
        return current_condition
    end
end

function parameter_condition_factory(config::Params.ConfigStruct)
    if config.pc == PC_LENGTH
        return length_condition
    elseif config.pc == PC_PERIOD
        return period_condition
    else
        throw(error("Unknown parameter criterion $(config.pc)"))
    end
end

end