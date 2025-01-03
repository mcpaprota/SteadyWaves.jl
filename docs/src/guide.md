# Package Guide

SteadyWaves is a small package for calculating parameters of nonlinear, steady water waves. In this guide, we present how to install and use the package. 

## Installation

```julia
pkg> add SteadyWaves
julia> using SteadyWaves
```

## Quick start

In this example, we set define basic parameters of a regular wave with respect to wavelength $L$ and wave height $H$, while we use [`CairoMakie.jl`](https://github.com/MakieOrg/Makie.jl) for graphical presentation of results
```@example 1
using SteadyWaves
using CairoMakie # use the plotting package

d = 1.0 # water depth (m)
H = 0.2 # wave height (m)
L = 2.0 # wavelength (m)
k = 2π / L # wave number (rad/s)
nothing # hide
```
and calculate wave profile using [`fourier_approx`](@ref) function along `2N+1` points with a parameter flag `pc=1`. Wave period `T` is calculated using [`wave_period`](@ref) function.
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
and calculate wave profile using [`fourier_approx`](@ref) function along `2N+1` points with a parameter flag `pc=2`. Wavelength is calculated using [`wavelength`](@ref) function.
```@example 2

N = 40 # set the number of points
u = fourier_approx(d, H, T; pc=2, cc=1, N=N) # apply Fourier Approximation Method
L = wavelength(u, d, N) # calculate wavelength (rad/s)
k = 2π / L # get wave number (rad/m)
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

We may also check the [`wave_height`](@ref).

```@example 2
println("The wave height is $(wave_height(u, d, N)) m.")
```

## Shoaling waves

Here, we present a more complicated example, which reproduces a shoaling diagram. We follow [Eldrup2020](@citet), who presented shoaling coefficient $K$ calculated by Fourier Approximation Method [Rienecker1981](@cite) as a function of water depth $d$ normalized by deep-water wave number $k_0$ (cf. Figure 2 in [Eldrup2020](@cite)).

First, we define test cases in terms of deep-water wave steepness $H_0/L_0$ and minimum values of water depth $d_{min}$ for each case.

```@example 3
using SteadyWaves
using CairoMakie

H₀_L₀ = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05] # set wave steepness cases
d_min = [0.019, 0.0275, 0.046, 0.0784, 0.1091, 0.1714] # set minimal depth values
nothing # hide
```

We set a number of points `N` for [`shoaling_approx`](@ref) (which uses in-place [`SteadyWaves.fourier_approx!`](@ref)), number of depth steps `N_d` and create a matrix container `results` for storing the outcome of analysis.

```@example 3
N = 40 # set number of points
N_d = 2000 # set number of steps for changing depth
results = zeros(N_d, 2 * length(H₀_L₀))
nothing # hide
```

We define initial deep-water conditions, where we calculate angular wave frequency `ω₀` according to linear dispersion relation.

```@example 3
# initial conditions
g = 9.81 # gravitational acceleration (m/s²)
d₀ = 1 # water depth (m)
L₀ = 2 # wavelength (m)
k₀ = 2π / L₀ # wave number (rad/m)
ω₀ = √(g * k₀ * tanh(k₀ * d₀)) # angular wave frequency (rad/s)
H₀ = H₀_L₀ * L₀ # wave heights (m)
nothing # hide
```

We loop over wave cases each time defining a vector of depth values `d` from deep to shallow water according to a predefined minimal value `d_min` and applying [`shoaling_approx`](@ref) function to calculate shoaling coefficient `K`. We store the results in `results`.

```@example 3
for i in eachindex(H₀)
    d = reverse(range(d_min[i], 1, N_d))  * d₀ # water depths (m)
    K = shoaling_approx(d, H₀[i], L₀; N=N) # shoaling coefficients
    # save results
    results[:, 2i-1] = k₀ * d
    results[:, 2i] = K
end
nothing # hide
```

For reference, we calculate linear shoaling coefficient

$$K = \sqrt{\frac{k\left(1+\frac{2k_0d_0}{\sinh2k_0d_0}\right)}{k_0\left(1+\frac{2kd}{\sinh2kd}\right)}}$$

which is independent of wave steepness $H_0/L_0$.

```@example 3
d = reverse(range(0.01, 1, N_d))  * d₀ # water depths (m)
k = wave_number.(d, ω₀) # linear wave numbers (rad/m)
# calculate linear shoaling coefficient
K = @. sqrt((k * (1 + 2k₀ * d₀ / sinh(2k₀ * d₀))) / (k₀ * (1 + 2k * d / sinh(2k * d))))
nothing # hide
```

Finally, we diagram the results.

```@example 3
# CairoMakie figure
with_theme(theme_latexfonts()) do # use latex theme
    fig = Figure(size = (400, 300), fontsize = 9)
    ax = Axis(fig[1, 1], 
        xlabel=L"$k_0d$",
        xscale=log10,
        xticks=[0.1, 1, 10],
        ylabel=L"$K$",
        yscale=log10,
        yticks=1:0.2:2.4,
        limits=(0.04, π, 0.9, 2.4);
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(10),
        xminorgridvisible=true,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(4),
        yminorgridvisible=true,
    )
    lines!(ax, k₀ * d, K;
        color=:tomato,
        linestyle=:dot,
        label="linear theory")
    for i in eachindex(H₀_L₀)
        lines!(ax, results[:, 2i-1], results[:, 2i];
            color=:dodgerblue4,
            linestyle=:dash, label="Rienecker and Fenton")
    end
    text!(ax, k₀ * d_min, [2.3, 2.0, 1.6, 1.35, 1.25, 1.1],
        text=[L"H_0/L_0= %$(value)" for value in H₀_L₀],
        offset=(5, 0),
        align=(:left, :center)
    )
    axislegend(ax, position=:rt, unique=true, patchsize=(42, 1))
    fig
end
```
