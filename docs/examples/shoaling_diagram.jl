using SteadyWaves
using CairoMakie

H₀toL₀

N_d = 2000 # set number of steps for changing depth
H_over_L = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05] # set wave stepness cases
# set minimal values for depth for wave stepness cases H_over_L
min_d₀ = [0.019, 0.0275, 0.046, 0.0784, 0.1091, 0.1714]
# create results matrix
results = zeros(n_steps_d, 2 * length(H_over_L))

# initial conditions
g = 9.81 # gravitational acceleration (m/s²)
d₀ = 1 # water depth (m)
L₀ = 2 # wavelength (m)
k₀ = 2π / L₀ # wave number (rad/m)
ω₀ = √(g * k₀ * tanh(k₀ * d₀)) # angular wave frequency (rad/s)
H = H_over_d * L₀ # wave height (m)
N = 40 # number of eigenvalues of solution

for i in eachindex(H)
    # initial nonlinear solution
    range_d = reverse(range(min_d₀[i], 1, N_d)) # water depths range (m)
    ds = range_d * d₀ # water depths (m)
    K = shoaling_approx(ds, H[i], L₀; N=N)
    results[:, 2i-1] = k₀ * ds
    results[:, 2i] = K
    println(i)
end

ds = reverse(range(0.01, 1, n_steps_d))
k = dispersion_relation.(ds, ω₀)

# initialize figure of size
fig_width = 5 * 72 # inches x points
aspect = 1.3
with_theme(theme_latexfonts()) do
    fig = Figure(size=(fig_width, fig_width / aspect), fontsize=9, figure_padding=4)
    ax = Axis(fig[1, 1],
        xlabel=L"$k_0d$",
        xscale=log10,
        xticks=[0.1, 1, 10],
        ylabel=L"$K_s$",
        yscale=log10,
        yticks=1:0.2:3.4,
        limits=(0.04, π, 0.9, 2.4);
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(10),
        xminorgridvisible=true,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(4),
        yminorgridvisible=true,
    )
    lines!(k₀ * ds, @. sqrt((k * (1 + 2k₀ * d₀ / sinh(2k₀ * d₀))) / (k₀ * (1 + 2k * ds / sinh(2k * ds))));
        color=:tomato,
        linestyle=:dot,
        label="linear theory")
    for n in 1:length(H_over_L)
        lines!(results[:, 2n-1], results[:, 2n];
            color=:dodgerblue4,
            linestyle=:dash, label="Rienecker and Fenton")
    end
    text!(ax, min_d₀ * k₀, [2.3, 2.0, 1.6, 1.35, 1.25, 1.1],
        text=[L"H_0/L_0= %$(H_L)" for H_L in H_over_L],
        offset=(5, 0),
        align=(:left, :center)
    )
    axislegend(ax, position=:rt, unique=true, patchsize=(42, 1))
    fig
    # save("figures/shoaling_diagram.png", fig, px_per_unit=10, pt_per_unit=1)
end
