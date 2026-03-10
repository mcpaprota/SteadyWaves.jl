# SPDX-License-Identifier: MIT

# Main module

"""
SteadyWaves is an implementation of Rienecker and Fenton (1981)
Fourier Approximation Method to steady, periodic, nonlinear waves
propagating in water of constant depth.
"""
module SteadyWaves

include("physics.jl")
using .Physics: G

include("params.jl")
using .Params

include("index.jl")

include("output.jl")
using .Output
using .Output: wave_height, wave_length, wave_power, wave_period,
    dimensionless_vertical_velocity, vertical_velocity,
    dimensionless_horizontal_velocity, horizontal_velocity

include("nonlinear_system.jl")

include("steady.jl")
using .Steady: fourier_approx, linear_wave_number

include("shoaling.jl")
using .Shoaling: topo_approx, fourier_approx!

export fourier_approx,fourier_approx!, topo_approx, wave_number
export wave_period, wave_length, wave_height, linear_wave_number
export dimensionless_vertical_velocity, vertical_velocity
export dimensionless_horizontal_velocity, horizontal_velocity
export CurrentCriterion, CC_EULER, CC_STOKES
export ParameterCriterion, PC_LENGTH, PC_PERIOD
export Output
export G

end
