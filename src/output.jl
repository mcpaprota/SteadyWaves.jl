# SPDX-License-Identifier: MIT

# Functions for calculating the output
module Output
using ..Physics

using ..Index

"""
    wave_period(u, d, N; g=G)

Calculate dimensional wave period `T` from solution `u`.
"""
function wave_period(u, d, N; g=G,id=Index.INDEX_STRUCT)
    T = 2π / u[id.C(N)] / √(u[id.D(N)]/ d * g)
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
    celerity = u[2N+C_INDEX]
    depth = u[2N+D_INDEX]
    flux = u[2N+Q_INDEX]
    velocity = u[2N+U_INDEX]
    bernoulli = u[2N+R_INDEX]

    relative_elevation = u[elevation_indexes(N)] .- depth

    U_e = celerity - velocity
    I_p = flux + depth * U_e
    E_p = (relative_elevation[1]^2 + relative_elevation[end]^2 + 2 * sum(relative_elevation[2:end-1] .^ 2)) / 4N
    Q = velocity / √depth - flux / depth^(1.5)
    E_k = 0.5 * (celerity * I_p - U_e * Q * depth^(1.5))
    U_b2 = 2 * bernoulli - celerity^2
    F = celerity * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + celerity * depth) + celerity * U_e * (velocity * depth - flux)
    return F
end

function stream_eigenfunction(hiperbolic,trigonometric,u,N,kx,kz,j)
    # takes index of jth stream function coefficient
    B_index = stream_indexes(N)[j] 

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

export C_INDEX, D_INDEX, H_INDEX, Q_INDEX, R_INDEX, U_INDEX

export elevation_indexes,stream_indexes

end