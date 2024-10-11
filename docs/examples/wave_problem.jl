# Schematics of a wave problem

using CairoMakie
using SteadyWaves

d = 1.0 # water depth (m)
H = 0.4 # wave height (m)
L = 10.0 # wavelength (m)
k = 2π / L # wave number (rad/s)

N = 40
u = fourier_approx(d, H, L; pc=1, cc=1, N=N)
x = range(0, L, 2N + 1)
η = [reverse(u[2:N+1]); u[1:N+1]]

# plotting
# initialize figure of size
fig_width = 5 * 72 # inches x points
aspect = 1.3
with_theme(theme_latexfonts()) do
    fig = Figure(size=(fig_width, fig_width / aspect), fontsize=9, figure_padding=4)
    ax = Axis(fig[1, 1],
        xlabel=L"$x/L$",
        ylabel=L"$k\eta$",
        limits=(0, L, -0.01, 1);
    )
lines!(ax, x, η) # wave profile
lines!(ax, x, x*0, color=:black) # horizontal bottom
lines!(ax, x, x*0 .+ k*d, color=:black, linestyle=:dash) # mean depth
lines!(ax, x, x*0 .+ η[1], color=:black, linewidth=0.5) # trough level
lines!(ax, x, x*0 .+ η[N+1], color=:black, linewidth=0.5) # crest level


fig
end
