# SPDX-License-Identifier: MIT

# Main module

"""
SteadyWaves is an implementation of Rienecker and Fenton (1981)
Fourier Approximation Method to steady, periodic, nonlinear waves
propagating in water of constant depth.
"""
module SteadyWaves

using NonlinearSolve

export fourier_approx, shoaling_approx

include("fam.jl")
include("shoaling.jl")

end
