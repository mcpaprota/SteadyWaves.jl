```@meta
CurrentModule = SteadyWaves
```

# SteadyWaves.jl

A Fourier Approximation Method to steady, periodic, nonlinear waves propagating in water of constant depth (Rienecker and Fenton, 1981).


## Overview

`SteadyWaves.jl` is a Julia package for calculation of properties of steady, periodic, and nonlinear waves within the framework of potential flow. The solution is derived using a Fourier Approximation Method applied to a periodic boundary value problem of waves propagating in water of constant depth up to a limiting height with additional capability to describe nonlinear shoaling waves (Rienecker and Fenton, 1981; Fenton, 1988, 1999). Unlike previous implementations of the method, here, we use [`NonlinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl) package to `solve()` a set of nonlinear equations using `RobustMultiNewton` mode.

## Wave problem

m
