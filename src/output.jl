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

"""
    wave_power(u, N)

Calculate dimensionless wave power `F` from solution `u` (non-public function).
"""
function wave_power(u, N)
    c = u[2N+C_INDEX]
    d = u[2N+D_INDEX]
    q = u[2N+Q_INDEX]
    v = u[2N+U_INDEX]
    r = u[2N+R_INDEX]

    relative_eta = u[eta_indexes(N)] .- d

    u_e = c - v
    I_p = q + d * u_e # mean wave momentum
    E_p = (relative_eta[1]^2 + relative_eta[end]^2 + 2 * sum(relative_eta[2:end-1] .^ 2)) / 4N # mean potential energy
    Q = v * d - q # volume flux
    E_k = 0.5 * (c * I_p - u_e * Q) # mean kinetic energy
    U_b2 = 2 * r - c^2 # mean square of bed velocity
    F = c * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + c * d) + c * u_e * Q # mean energy flux - wave power
    return F
end

function stream_eigenfunction(hiperbolic,trigonometric,u,N,kx,kz,j)
    # takes index of jth stream function coefficient
    B_index = psi_indexes(N)[j] 

    # retrieves jth stream function coefficient without creating an array
    B = u[B_index]

    return B * hiperbolic(j * kz) / cosh(j * u[2N+D_INDEX]) * trigonometric(j * kx)
end

function surface_stream_eigenfunction(hiperbolic,trigonometric,u,N,m,j)
   return stream_eigenfunction(hiperbolic, trigonometric, u, N, m/N * π, u[m+1], j) 
end

function dimensionless_vertical_velocity(u,N,kx,kz)
    return sum([j*stream_eigenfunction(sinh, sin, u, N, kx, kz, j) for j in 1:N])
end


function vertical_velocity(u,N,x,z,k,g=G)
    return sqrt(g/k)*dimensionless_vertical_velocity(u,N,k *x, k * z)
end

function dimensionless_horizontal_velocity(u,N,kx,kz)
    return -u[2N + U_INDEX] + sum([j*stream_eigenfunction(cosh, cos, u, N, kx, kz, j) for j in 1:N])
end


function horizontal_velocity(u,N,x,z,k,g=G)
    return sqrt(g/k)*dimensionless_horizontal_velocity(u,N,k *x, k * z)
end

function dimensionless_pressure(u,N,kx,kz)
    v_x = dimensionless_horizontal_velocity(u,N,kx,kz)
    v_z = dimensionless_vertical_velocity(u,N,kx,kz)

    return u[2N+R_INDEX] - v_x^2 / 2 - v_z^2 / 2 - kz + u[2N+D_INDEX]
end

function pressure(u,N,x,z,k,g=9.81,rho=RHO)
    return rho * g / k * dimensionless_pressure(u, N, k * x, k * z)
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
