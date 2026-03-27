# SPDX-License-Identifier: MIT

# Functions for calculating the output
module Output
using ..Physics

using ..Index
using ..Wave: WaveStruct

"""
    wave_period(u, d, N; g=G)

Calculate dimensional wave period `T` from solution `u`.
"""
function wave_period(u, d, N; g=G)
    T = 2π / u[2N+C_INDEX] / √(u[2N+D_INDEX] / d * g)
    return T
end

"""
    wavelength(u, d, N)

Calculate dimensional wavelength `L` from solution `u`.
"""
function wavelength(u, d, N)
    L = 2π * d / u[2N+D_INDEX]
    return L
end

"""
    wave_number(u, d, N)

Calculate dimensional wave number `k` from solution `u`.
"""
function wave_number(u, d, N)
    return u[2N+D_INDEX] / d
end

"""
    wave_height(u, d, N)

Calculate dimensional wave height `H` from solution `u`.
"""
function wave_height(u, d, N)
    H = u[2N+H_INDEX] / u[2N+D_INDEX] * d
    return H
end

function vertical_velocity(w,kx,kz)
    return w.v.z(kx,kz)
end

function vertical_velocity(w, x, z, df)
    return w.v.z(x * df.L,z * df.D) / df.v.z
end

function horizontal_velocity(w,kx,kz)
    return w.v.x(kx,kz)
end

function horizontal_velocity(w, x, z, df)
    return w.v.x(x * df.L,z * df.D) / df.v.x
end

function indirect_pressure(w,kx,kz)
    return w.R - w.v.x(kx,kz)^2 / 2 - w.v.z(kx,kz)^2 / 2 - kz + w.D
end

function pressure(w,kx,kz)
    return w.P === nothing ? indirect_pressure(w,kx,kz) : w.P(kx,kz)   
end

function pressure(w,x,z,df)
    return pressure(w,x * df.L, z * df.D) / df.P    
end

function indirect_wave_period(w::WaveStruct)
    return w.L / w.C
end

function wave_period(w::WaveStruct)
    return w.T === nothing ? indirect_wave_period(w) : w.T
end

function wave_period(w::WaveStruct, df::WaveStruct)
    return wave_period(w) / df.T    
end

function indirect_wavelength(w::WaveStruct)
    return w.T * w.C
end

function wavelength(w::WaveStruct)
    return w.L === nothing ? indirect_wavelength(w) : w.L
end

function wavelength(w::WaveStruct,df::WaveStruct)
    return wavelength(w) / df.L
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

function wave_power(w::WaveStruct)
    return w.F === nothing ? indirect_wave_power(w) : w.F
end

function wave_power(w::WaveStruct,df::WaveStruct)
    return wave_power(w) / df.F
end

export C_INDEX, D_INDEX, H_INDEX, Q_INDEX, R_INDEX, U_INDEX

export eta_indexes,psi_indexes

end
