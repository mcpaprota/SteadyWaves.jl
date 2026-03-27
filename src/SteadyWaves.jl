# SPDX-License-Identifier: MIT

# Main module

"""
SteadyWaves is an implementation of Rienecker and Fenton (1981)
Fourier Approximation Method to steady, periodic, nonlinear waves
propagating in water of constant depth.
"""
module SteadyWaves

include("physics.jl")
using .Physics: G, RHO

include("params.jl")
using .Params

include("index.jl")

include("surface.jl")

include("wave.jl")

include("linear.jl")
using .Linear: linear_wave_number, linear_solution

include("output.jl")
using .Output
using .Output: wave_height, wavelength, wave_power, wave_period,
    dimensionless_vertical_velocity, vertical_velocity,
    dimensionless_horizontal_velocity, horizontal_velocity,
    dimensionless_pressure, pressure

include("condition.jl")

include("nonlinear_system.jl")

include("steady.jl")
using .Steady: fourier_approx

include("shoaling.jl")
using .Shoaling: topo_approx, update_depth_fourier_approx

export fourier_approx,update_depth_fourier_approx, topo_approx, wave_number
export wave_period, wavelength, wave_height, linear_wave_number
export dimensionless_vertical_velocity, vertical_velocity
export dimensionless_horizontal_velocity, horizontal_velocity
export dimensionless_pressure, pressure
export CurrentCriterion, CC_EULER, CC_STOKES
export ParameterCriterion, PC_LENGTH, PC_PERIOD
export Output
export G, RHO

end
