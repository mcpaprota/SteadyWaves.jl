# SPDX-License-Identifier: MIT

# Functions for calculating the output
module Output

using ..Index

"""
    wave_period(u, d, N; g=9.81)

Calculate dimensional wave period `T` from solution `u`.
"""
function wave_period(u, d, N; g=9.81)
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

export C_INDEX, D_INDEX, H_INDEX, Q_INDEX, R_INDEX, U_INDEX

export elevation_indexes,stream_indexes

end