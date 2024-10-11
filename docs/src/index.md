```@meta
CurrentModule = SteadyWaves
```

# SteadyWaves.jl

A Fourier Approximation Method to steady, periodic, nonlinear waves propagating in water of constant depth [Rienecker1981](@cite).


## Overview

`SteadyWaves.jl` is a Julia package for calculation of properties of steady, periodic, and nonlinear waves within the framework of potential flow. The solution is derived using a Fourier Approximation Method applied to a periodic boundary value problem of waves propagating in water of constant depth up to a limiting height [Rienecker1981,Fenton1988,Fenton1999](@cite) with additional capability to describe nonlinear shoaling waves [Rienecker1981](@cite). Unlike previous implementations of the method, here, we use [`NonlinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl) package [Pal2024](@cite) to `solve` a set of nonlinear equations using `RobustMultiNewton` mode.

## Wave problem

We consider steady, periodic waves propagating in water of constant depth. In figure, a schematic view of the wave problem is presented with coordinate system definition.

```@example
include("../examples/wave_problem.jl") # hide
```