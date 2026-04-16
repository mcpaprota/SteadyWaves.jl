# SPDX-License-Identifier: MIT

# Main module

"""
SteadyWaves is an implementation of Rienecker and Fenton (1981)
Fourier Approximation Method to steady, periodic, nonlinear waves
propagating in water of constant depth.
"""
module SteadyWaves

include("wave/physics.jl")
using .Physics: G, RHO

include("system/params.jl")
using .Params

include("wave/index.jl")

include("wave/surface.jl")

include("wave/velocity.jl")

include("wave/wave.jl")

include("crapper.jl")

include("wave/dimensional_factor.jl")

include("linear.jl")
using .Linear: linear_wave_number, linear_solution

include("output.jl")
using .Output
using .Output: wave_height, wavelength, wave_power, wave_period,
    vertical_velocity, horizontal_velocity, pressure, elevation

include("condition.jl")

include("system/nonlinear_system.jl")

include("steady.jl")
using .Steady: fourier_approx

include("shoaling.jl")
using .Shoaling: topo_approx, update_depth_fourier_approx



export fourier_approx,update_depth_fourier_approx, topo_approx, wave_number
export wave_period, wavelength, wave_height, linear_wave_number, elevation
export vertical_velocity
export horizontal_velocity
export pressure
export CurrentCriterion, CC_EULER, CC_STOKES
export ParameterCriterion, PC_LENGTH, PC_PERIOD

export Output
export Params
export Physics
export NonlinearSystem
export Steady
export Shoaling
export Index
export Linear
export Crapper

export G, RHO

end
