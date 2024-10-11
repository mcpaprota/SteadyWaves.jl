# Package Guide

SteadyWaves is a small package for calculating parameters of nonlinear, steady water waves. In this guide, we present how to install and use the package. 

## Installation

```julia
pkg> add SteadyWaves
julia> using SteadyWaves
```

## Quick start

```@setup 1
using SteadyWaves
```

In this example, we set basic parameters of a regular wave
```@example 1
d = 1.0 # water depth (m)
H = 0.2 # wave height (m)
L = 2.0 # wavelength (m)
k = 2π / L # wave number (rad/s)
nothing # hide
```
and calculate its profile using [`fourier_approx`](@ref) function along `2N+1` points
```@example 1
# use the plotting package
using CairoMakie

N = 40 # set the number of points
u = fourier_approx(d, H, L; pc=1, cc=1, N=N) # apply Fourier Approximation Method
kη = [reverse(u[2:N+1]); u[1:N+1]] # vcat non-dimensional profile vector and its reverse
x = range(0, L, 2N + 1) # discretize L to match kη
η = kη / k # use dimensional values

# CairoMakie figure
with_theme(theme_latexfonts()) do # use latex theme
    fig = Figure(size = (300, 200), fontsize = 9)
    ax = Axis(fig[1,1], 
            xlabel=L"$x$ (m)",
            ylabel=L"$\eta$ (m)",
            limits=(0, L, 0, η[N+1] * 1.1);) # set upper limit to 110% of η_max
    lines!(ax, x, η)
    fig
end
```