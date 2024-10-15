# Schematics of a wave problem

using CairoMakie
using SteadyWaves

include("color_scheme.jl")

# setting the wave

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
aspect = 1.5
with_theme(theme_latexfonts()) do
    fig = Figure(size=(fig_width, fig_width / aspect), fontsize=9, figure_padding=4)
    ax = Axis(fig[1, 1],
        xlabel=L"$x/L$",
        ylabel=L"$k\eta$",
        limits=(0, L, -0.05, 0.85);
    )
    lines!(ax, x, η, color=water_surf) # wave profile
    band!(ax, x, η, 0, color=water_bulk)
    lines!(ax, x, x * 0, color=sand_surf) # horizontal bottom
    band!(ax, x, x * 0, -0.05, color=sand_bulk)
    lines!(ax, x, x * 0 .+ k * d, color=:black, linestyle=:dash, linewidth=0.9) # mean depth
    lines!(ax, x[1:30], x[1:30] * 0 .+ η[1], color=:black, linewidth=0.5) # trough level
    lines!(ax, x[20:N+1], x[20:N+1] * 0 .+ η[N+1], color=:black, linewidth=0.5) # crest level
    lines!(ax, [0.7L, 0.8L], [0.8, 0.8], color=:black, linewidth=1)

    # coordinate systems and arrows
    coord_sys_shift = 0.4L
    n_eta = 50 # eta symbol position
    lines!(ax, [L / 2, 0.68L] .- coord_sys_shift, [0, 0], color=:black, linewidth=0.5) # x axis
    lines!(ax, [0.5L, 0.5L] .- coord_sys_shift, [0, 0.25d], color=:black, linewidth=0.5) # z axis
    points_axis = [(L / 2 - coord_sys_shift, 0.25 * d), (0.68L - coord_sys_shift, 0)]
    scatter!(ax, points_axis; color=:black,
        marker=:utriangle, markersize=(6, 12), rotation=[0, -π / 2]) # arrows
    scatter!(ax, (0.5L - coord_sys_shift, 0), strokewidth=0.5, markersize=6, color=:gray)
    lines!(ax, [0.3L, 0.3L], [η[1], η[N+1]], color=:black, linewidth=0.5) # height
    lines!(ax, [0.8L, 0.8L], [0, k * d], color=:black, linewidth=0.5) # depth
    lines!(ax, [0, L], [0.44d, 0.44d], color=:black, linewidth=0.5) # length
    m_shift = 0.022
    points_arrows = [(0.3L, η[1] + m_shift), (0.3L, η[N+1] - m_shift),
        (0.8L, m_shift), (0.8L, k * d - m_shift), (8m_shift, 0.44d),
        (L - 8m_shift, 0.44d), (0.8L, 0.8)]
    scatter!(ax, points_arrows; color=:black,
        marker=:utriangle, markersize=(6, 12), rotation=[-π, 0, -π, 0, π / 2, -π / 2, -π / 2]) # arrows
    points_text = [(L / 1.95 - coord_sys_shift, 0.26 * d), (0.7L - coord_sys_shift, 0),
        (0.28L, (η[1] + η[N+1]) / 2), (0.78L, k * d / 2), (x[n_eta], η[n_eta]),
        (L / 2, 0.45d), (0.75L, 0.81)]
    text!(ax, points_text;
        text=[L"$z$", L"$x$", L"$H$", L"$d$", L"$\eta(x, t)$", L"$L$", L"$c$"],
        align=[(:left, :bottom), (:left, :top), (:right, :bottom),
            (:right, :bottom), (:left, :bottom), (:center, :bottom), (:center, :bottom)])
    hidedecorations!(ax)
    hidespines!(ax)
    save("docs/src/figures/wave_problem.png", fig, px_per_unit=10, pt_per_unit=1)
end
