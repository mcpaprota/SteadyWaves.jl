module Indirect

import ..Wave: WaveStruct
import ..Params: Params, ConfigStruct

function indirect_wave_period(w)
    return w.L / w.C
end

function indirect_wavelength(w)
    return w.C * w.T
end

function indirect_dynamic_pressure(w,kx,kz)
    return w.R - w.v.x(kx,kz)^2 / 2 - w.v.z(kx,kz)^2 / 2
end

function indirect_pressure(w,kx,kz)
    return w.R - w.v.x(kx,kz)^2 / 2 - w.v.z(kx,kz)^2 / 2 - kz + w.D
end

function indirect_pressure(w,kx,kz,u)
    return w.R - (w.v.x(kx,kz)+u)^2 / 2 - w.v.z(kx,kz)^2 / 2 - kz + w.D
end

function indirect_wave_period(w::WaveStruct)
    return w.L / w.C
end

function indirect_wavelength(w::WaveStruct)
    return w.T * w.C
end

function indirect_celerity_from_euler_condition(w::WaveStruct)
    return w.U
end

function indirect_celerity_from_stokes_condition(w::WaveStruct)
    return w.U - w.Q / w.D
end

function indirect_celerity_factory(config::ConfigStruct)
    if config.cc == Params.CC_STOKES
        return indirect_celerity_from_stokes_condition
    elseif config.cc == Params.CC_EULER
        return indirect_celerity_from_euler_condition
    else
        throw(error("Unknown current criterion $(config.cc)"))
    end
end

function indirect_wave_power(w)
    u_e = w.C - w.U
    I_p = w.Q + w.D * u_e # mean wave momentum
    Q = w.U * w.D - w.Q # volume flux
    E_k = 0.5 * (w.C * I_p - u_e * Q) # mean kinetic energy
    U_b2 = 2 * w.R - w.C^2 # mean square of bed velocity
    F = w.C * (3E_k - 2w.eta.e_p) + 0.5 * U_b2 * (I_p + w.C * w.D) + w.C * u_e * Q # mean energy flux - wave power
    return F
end

function indirect_surface_tension(sigma,kz_d1,kz_d2)
    return sigma * kz_d2 / (1 + kz_d1^2)^1.5
end

end