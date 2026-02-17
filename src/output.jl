# SPDX-License-Identifier: MIT

# Functions for calculating the output

"""
    wave_period(u, d, N; g=9.81)

Calculate dimensional wave period `T` from solution `u`.
"""
function wave_period(u, d, N; g=9.81)
    T = 2π / u[2N+2] / √(u[2N+3] / d * g)
    return T
end

"""
    wavelength(u, d, N)

Calculate dimensional wavelength `L` from solution `u`.
"""
function wavelength(u, d, N)
    L = 2π * d / u[2N+3]
    return L
end

"""
    wave_height(u, d, N)

Calculate dimensional wave height `H` from solution `u`.
"""
function wave_height(u, d, N)
    H = u[2N+7] / u[2N+3] * d
    return H
end

"""
    wave_power(u, N)

Calculate dimensionless wave power `F` from solution `u` (non-public function).
"""
function wave_power(u, N)
    U_e = u[2N+2] - u[2N+6]
    I_p = u[2N+4] + u[2N+3] * U_e
    E_p = ((u[1] - u[2N+3])^2 + (u[N+1] - u[2N+3])^2 + 2 * sum((u[2:N] .- u[2N+3]) .^ 2)) / 4N
    Q = u[2N+6] / √u[2N+3] - u[2N+4] / u[2N+3]^(1.5)
    E_k = 0.5 * (u[2N+2] * I_p - U_e * Q * u[2N+3]^(1.5))
    U_b2 = 2u[2N+5] - u[2N+2]^2
    F = u[2N+2] * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + u[2N+2] * u[2N+3]) + u[2N+2] * U_e * (u[2N+6] * u[2N+3] - u[2N+4])
    return F
end

"""
    wave_number(d, ω, g=9.81, ϵ=10^-12)

Calculate wave_number `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""
function wave_number(d, ω, g=9.81, ϵ=10^-12)
    k = k₀ = ω^2 / g # initial guess
    while max(abs(k * tanh(k * d) - k₀)) > ϵ
        k = k₀ / tanh(k * d)
    end
    return k
end