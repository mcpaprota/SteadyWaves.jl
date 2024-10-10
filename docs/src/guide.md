# Package Guide

SteadyWaves is a small package for calculating parameters of nonlinear, steady water waves. In this guide, we present how to install and use the package. 

## Installation

```julia
pkg> add SteadyWaves
julia> using SteadyWaves
```

## Quick start

```@setup 1
#using Pkg
#Pkg.add("CairoMakie")
#Pkg.develop(url="https://github.com/mcpaprota/SteadyWaves.jl")
#using CairoMakie
using CairoMakie
using SteadyWaves
```

In this example, we set basic parameters of a regular wave
```@example 1
d = 1.0 # water depth (m)
H = 0.1 # wave height (m)
L = 1.0 # wavelength (m)
nothing # hide
```
and calculate its half-profile using [`fourier_approx`](@ref) function along `N` points
```@example 1
N = 20
u = fourier_approx(d, H, L; pc=1, cc=1, N=N)
η = u[1:N+1]
lines(η)
```