using SteadyWaves
using CairoMakie

H₀_L₀ = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05] # set wave steepness cases
d_min = [0.019, 0.0275, 0.046, 0.0784, 0.1091, 0.1714] # set minimal depth values

N = 40 # set number of points
N_d = 600 # set number of steps for changing depth
results = zeros(N_d, 2 * length(H₀_L₀))

# initial conditions
g = 9.81 # gravitational acceleration (m/s²)
d₀ = 1 # water depth (m)
L₀ = 2 # wavelength (m)
k₀ = 2π / L₀ # wave number (rad/m)
ω₀ = √(g * k₀ * tanh(k₀ * d₀)) # angular wave frequency (rad/s)
H₀ = H₀_L₀ * L₀ # wave heights (m)

for i in eachindex(H₀)
    d = reverse(logrange(d_min[i], 1, N_d))  * d₀ # water depths (m)
    K = topo_approx(d, H₀[i], L₀; N=N) # shoaling coefficients
    # save results
    results[:, 2i-1] = k₀ * d
    results[:, 2i] = K
end

d = reverse(logrange(0.01, 1, N_d))  * d₀ # water depths (m)
k = wave_number.(d, ω₀) # linear wave numbers (rad/m)
# calculate linear shoaling coefficient
K = @. sqrt((k * (1 + 2k₀ * d₀ / sinh(2k₀ * d₀))) / (k₀ * (1 + 2k * d / sinh(2k * d))))

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
