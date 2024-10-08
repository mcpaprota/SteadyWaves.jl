"""
SteadyWaves is an implementation of Rienecker and Fenton (1981)
Fourier Approximation Method to steady, periodic, nonlinear waves
propagating in water of constant depth.
"""
module SteadyWaves

using NonlinearSolve

export fourier_approx

include("functions.jl")

end
