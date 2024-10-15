# Package Guide

SteadyWaves is a small package for calculating parameters of nonlinear, steady water waves. In this guide, we present how to install and use the package. 

## Installation

```julia
pkg> add SteadyWaves
julia> using SteadyWaves
```

## Quick start

In this example, we set define basic parameters of a regular wave with respect to wavelength $L$ and wave height #H#, while we use CairoMakie for graphical presentation of results
```@example 1
using SteadyWaves
using CairoMakie # use the plotting package

d = 1.0 # water depth (m)
H = 0.2 # wave height (m)
L = 2.0 # wavelength (m)
k = 2π / L # wave number (rad/s)
nothing # hide
```
and calculate wave profile using [`fourier_approx`](@ref) function along `2N+1` points with a parameter flag `pc=1`
```@example 1

N = 40 # set the number of points
u = fourier_approx(d, H, L; pc=1, cc=1, N=N) # apply Fourier Approximation Method
kη = [reverse(u[2:N+1]); u[1:N+1]] # vcat non-dimensional profile vector and its reverse
T = wave_period(u, d, N) # calculate wave period
x = range(0, L, 2N + 1) # discretize L to match kη
η = kη / k # use dimensional values

# CairoMakie figure
with_theme(theme_latexfonts()) do # use latex theme
    fig = Figure(size = (300, 200), fontsize = 9)
    ax = Axis(fig[1,1], 
            xlabel=L"$x$ (m)",
            ylabel=L"$\eta$ (m)",
            # display wave parameters in the title
            title=L"$d=%$(d)$ m, $L=%$(L)$ m, $T=%$(round(T, digits=2))$ s, $H=%$(H)$ m",
            limits=(0, L, 0, η[N+1] * 1.1);) # set upper limit to 110% of η_max
    lines!(ax, x, η)
    fig # display figure
end
```

In another simple example, we set basic parameters of a regular wave with respect to wave period $T$ and height $H$
```@example 2
using SteadyWaves
using CairoMakie

d = 1.0 # water depth (m)
H = 0.2 # wave height (m)
T = 2.0 # wave period (s)
nothing # hide
```
and calculate wave profile using [`fourier_approx`](@ref) function along `2N+1` points with a parameter flag `pc=2`
```@example 2

N = 40 # set the number of points
u = fourier_approx(d, H, T; pc=2, cc=1, N=N) # apply Fourier Approximation Method
k = u[2N+3] / d # get wave number
L = 2π / k # calculate wavelength (rad/s)
kη = [reverse(u[2:N+1]); u[1:N+1]] # vcat non-dimensional profile vector and its reverse
x = range(0, L, 2N + 1) # discretize L to match kη
η = kη / k # use dimensional values

# CairoMakie figure
with_theme(theme_latexfonts()) do # use latex theme
    fig = Figure(size = (300, 200), fontsize = 9)
    ax = Axis(fig[1,1], 
            xlabel=L"$x$ (m)",
            ylabel=L"$\eta$ (m)",
            # display wave parameters in the title
            title=L"$d=%$(d)$ m, $L=%$(round(L, digits=2))$ m, $T=%$(T)$ s, $H=%$(H)$ m",
            limits=(0, L, 0, η[N+1] * 1.1);) # set upper limit to 110% of η_max
    lines!(ax, x, η)
    fig # display figure
end
```

## Shoaling waves

Here, we present a more complicated example, which reproduces a shoaling diagram. We follow [Eldrup2020](@citet), who presented shoaling coefficient $K_s$ calculated by Fourier Approximation Method [Rienecker1981](@cite) as a function water depth $d$ normalized by deep-water wave number $k_0$ (cf. Figure 2 in [Eldrup2020](@cite)).

First, we define test cases in terms of deep-water wave steepness $H_0/L_0$ and minimum values of water depth $d_{min}$ for each case.

```@example 2
using SteadyWaves
using CairoMakie

H₀_L₀ = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05] # set wave steepness cases
d_min = [0.019, 0.0275, 0.046, 0.0784, 0.1091, 0.1714] # set minimal depth values
nothing # hide
```

We set a number of points `N` for `fourier_approx`, number of depth steps `N_d` and create a matrix container `results` for storing the outcome of analysis.

```@example 2
N = 40 # set number of points
N_d = 2000 # set number of steps for changing depth
results = zeros(n_steps_d, 2 * length(H₀_L₀))
nothing # hide
```

We define initial deep-water conditions, where we calculate angular wave frequency `ω₀` according to linear dispersion relation.

```@example 2
# initial conditions
d₀ = 1 # water depth (m)
L₀ = 2 # wavelength (m)
k₀ = 2π / L₀ # wave number (rad/m)
ω₀ = √(g * k₀ * tanh(k₀ * d₀)) # angular wave frequency (rad/s)
H₀ = H₀_L₀ * L₀ # wave heights (m)
nothing # hide
```

We loop over wave cases each time defining a vector of depth values `d` from deep to shallow water according to a predefined minimal value `d_min` and applying `shoaling approx` function to calculate shoaling coefficient `K`. We store the results in `results`.

```@example 2
for i in eachindex(H)
    d = reverse(range(d_min[i], 1, N_d))  * d₀ # water depths (m)
    K = shoaling_approx(d, H₀[i], L₀; N=N) # 
    results[:, 2i-1] = k₀ * d
    results[:, 2i] = K
end
nothing # hide
```

For reference, we calculate linear shoaling coefficient (see e.g. eq.), which is independent of wave steepness $H_0/L_0$.

```@example 2
d = reverse(range(0.01, 1, N_d))  * d₀ # water depths (m)
k = dispersion_relation.(d, ω₀) # linear wave numbers (rad/m)
K = @. sqrt((k * (1 + 2k₀ * d₀ / sinh(2k₀ * d₀))) / (k₀ * (1 + 2k * d / sinh(2k * d)))) # linear shoaling coefficient
nothing # hide
```
