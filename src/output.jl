# SPDX-License-Identifier: MIT

# Functions for calculating the output

"""
`u[elevation_indexes(N)]`: free surface elevations *kη*
"""
function elevation_indexes(N)
    return 1:N+1
end

"""
`u[stream_indexes(N)]` : stream function coefficients *B*
"""
function stream_indexes(N)
    return N+2:2N+1
end

"""
`u[2N+C_INDEX]`: wave celerity *c√(k/g)*
"""
const C_INDEX::Int = 2

"""
`u[2N+D_INDEX]`: mean water depth *kη̄*
"""
const D_INDEX::Int = 3

"""
`u[2N+Q_INDEX]`: volume flux due to waves *q√(k³/g)*
"""
const Q_INDEX::Int = 4

"""
`u[2N+R_INDEX]`: Bernoulli constant *rk/g*
"""
const R_INDEX::Int = 5

"""
`u[2N+U_INDEX]`: mean flow velocity *Ū√(k/g)*
"""
const U_INDEX::Int = 6

"""
`u[2N+H_INDEX]`: wave height *kH*
"""
const H_INDEX::Int = 7


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

function reflect(array)
    return [reverse(array); array[2:end]]
end

function free_surface(u,N)
    return reflect(u[1:N+1])
end

function stream_eigenfunction(hiperbolic,trigonometric,u,N,kx,ky,j)
    return u[N+1+j] * hiperbolic(j * ky) / cosh(j * u[2N+3]) * trigonometric(j * kx * π / N)
end

function surface_stream_eigenfunction(hiperbolic,trigonometric,u,N,m,j)
   return stream_eigenfunction(hiperbolic,trigonometric,u,N,m,u[m+1],j) 
end


function u(u,N,kx,ky)
    k = 1 #TODO K wave_number
    return u[N+1] + k*sum([j*stream_eigenfunction(cosh,cos,u,N,kx,ky,j) for j in 1:N])
end