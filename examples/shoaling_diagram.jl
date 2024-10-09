n_steps_d = 5000
min_d₀ = [0.019, 0.0275, 0.046, 0.0784, 0.1091, 0.1714]
range_H_L = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05]
results = zeros(n_steps_d, 2 * length(range_H_L))

# initial conditions
d₀ = 1 # water depth (m)
L₀ = 2 # wavelength (m)
k₀ = 2π / L₀
ω₀ = √(g * k₀ * tanh(k₀ * d₀))
Hs = range_H_L * L₀ # wave height (m)
N = 40 # eigenvalues of solution

for i in eachindex(Hs)
    # initial nonlinear solution
    range_d = reverse(range(min_d₀[i], 1, n_steps_d)) # water depths (m)
    ds = range_d * d₀ # water depths (m)
    u = fam_L(d₀, Hs[i], L₀, 1, N, 1, false)
    H₀ = u[2N+7] / u[2N+3] * d₀
    F = u[2N+8]
    T = u[2N+9]
    println((F, T))
    # subsequent solutions at new depths
    for j in eachindex(ds)
        j == 1 ? fams2!(u, ds[1], 1, F, T, 1, N) : fams2!(ds[j], ds[j-1], F, T, 1, N, u)
        results[j, 2i-1] = k₀ * ds[j]
        results[j, 2i] = u[2N+7] / u[2N+3] * ds[j] / H₀
        println((i, j, ds[j], u[2N+7] / u[2N+3] * ds[j] / H₀, u[2N+8], u[2N+9]))
    end
end

ds = reverse(range(0.01, 1, n_steps_d))
k = disp_rel.(ds, ω₀)
