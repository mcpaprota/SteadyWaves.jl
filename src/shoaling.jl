# SPDX-License-Identifier: MIT

# Functions for shoaling calculations based on Fourier Approximation Method

"""
    fams!(kd, F, τ, N, M, u, fs_amps)

Get subsequent shoaling wave solution from Fourier Approximation Method
for mean water depths `d`, wave power `F`, wave period `τ`, number of eigenvalues `N`,
    and number of height steps `M`. Flags: cc - current criterion 1 - stokes, 2 - euler;
    fs_amps Fourier amplitudes or elevation

"""
function fourier_approx!(u, d, d_p, F, T; cc=1, N=10, g=9.81)
    init_conditions!(d_p / d, u, N)
    params = [F / (√g^3 * √d^5), T * sqrt(g / d), cc]
    problem = NonlinearProblem(f_1, u[1:2N+7], params)
    solution = solve(problem, RobustMultiNewton())
    u[1:2N+7] = sol.u
    U_e = u[2N+2] - u[2N+6]
    I_p = u[2N+4] + u[2N+3] * U_e
    E_p = ((u[1] - u[2N+3])^2 + (u[N+1] - u[2N+3])^2 + 2 * sum((u[2:N] .- u[2N+3]) .^ 2)) / 4N
    Q = u[2N+6] / √u[2N+3] - u[2N+4] / u[2N+3]^(1.5)
    E_k = 0.5 * (u[2N+2] * I_p - U_e * Q * u[2N+3]^(1.5))
    U_b2 = 2u[2N+5] - u[2N+2]^2
    F = u[2N+2] * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + u[2N+2] * u[2N+3]) + u[2N+2] * U_e * (u[2N+6] * u[2N+3] - u[2N+4])
    u[2N+8] = F * √(g^3/(u[2N+3]/d)^5) # add F√(k⁵/g³)/ρ to the outcome of the solution
    u[2N+9] = 2π / √(u[2N+3] / d * g) / u[2N+2] # add T to the outcome of the solution
    return u
end

function init_conditions!(ratio_d, u, N)
    u[1:N+1] =  1 .+ (u[1:N+1] .- 1) / ratio_d # kη
    u[N+2:2N+1] = u[N+2:2N+1] / √ratio_d # B
    u[2N+2] = u[2N+2] / √ratio_d # c√(k/g)
    u[2N+3] = u[2N+3] / ratio_d # kη̄
    u[2N+4] = u[2N+4] / √ratio_d^3 # q√(k³/g)
    u[2N+5] = u[2N+5] / ratio_d # rk/g
    u[2N+6] = u[2N+6] / √ratio_d # Ū√(k/g)
    u[2N+7] = u[2N+7] / ratio_d # kH
    return u
end

function f_1(du, u, p)
    N = (length(u) - 7) ÷ 2
    for m in 0:N
        Σ₁ = sum([u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N) for j in 1:N])
        du[m+1] = Σ₁ - u[2N+6] * (u[m+1] - u[2N+3]) - u[2N+4]
        Σ₂ = sum([j * u[N+1+j] * cosh(j * u[m+1]) / cosh(j * u[2N+3]) * cos(j * m * π / N) for j in 1:N])
        Σ₃ = sum([j * u[N+1+j] * sinh(j * u[m+1]) / cosh(j * u[2N+3]) * sin(j * m * π / N) for j in 1:N])
        du[N+1+m+1] = (-u[2N+6] + Σ₂)^2 / 2 + Σ₃^2 / 2 + u[m+1] - u[2N+3] - u[2N+5]
    end
    du[2N+3] = ((u[1] + u[N+1]) / 2 + sum(u[2:N])) / N - u[2N+3]
    du[2N+4] = u[1] - u[N+1] - u[2N+7]
    du[2N+5] = u[2N+2] * p[2] * √u[2N+3] - 2π
    p[3] == 1 ? du[2N+6] = u[2N+6] - u[2N+2] - u[2N+4] / u[2N+3] : du[2N+6] = u[2N+6] - u[2N+2]
    U_e = u[2N+2] - u[2N+6]
    I_p = u[2N+4] + u[2N+3] * U_e
    E_p = ((u[1] - u[2N+3])^2 + (u[N+1] - u[2N+3])^2 + 2 * sum((u[2:N] .- u[2N+3]) .^ 2)) / 4N
    Q = u[2N+6] / √u[2N+3] - u[2N+4] / u[2N+3]^(1.5)
    E_k = 0.5 * (u[2N+2] * I_p - U_e * Q * u[2N+3]^(1.5))
    U_b2 = 2u[2N+5] - u[2N+2]^2
    du[2N+7] = u[2N+2] * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + u[2N+2] * u[2N+3]) + u[2N+2] * U_e * (u[2N+6] * u[2N+3] - u[2N+4]) - p[1] * √u[2N+3]^5
    return nothing
end

function wave_power(u, N)
    U_e = u[2N+2] - u[2N+6]
    I_p = u[2N+4] + u[2N+3] * U_e
    E_p = ((u[1] - u[2N+3])^2 + (u[N+1] - u[2N+3])^2 + 2 * sum((u[2:N] .- u[2N+3]) .^ 2)) / 4N
    Q = u[2N+6] / √u[2N+3] - u[2N+4] / u[2N+3]^(1.5)
    E_k = 0.5 * (u[2N+2] * I_p - U_e * Q * u[2N+3]^(1.5))
    U_b2 = 2u[2N+5] - u[2N+2]^2
    F = u[2N+2] * (3E_k - 2E_p) + 0.5 * U_b2 * (I_p + u[2N+2] * u[2N+3]) + u[2N+2] * U_e * (u[2N+6] * u[2N+3] - u[2N+4])
    return F
end

F = wave_power(u, N) # calculate wave power F
push!(u, 2π / L * H) # add kH to the outcome of the solution
push!(u, F * √(g^3/(2π/L)^5)) # add F/ρ to the outcome of the solution
push!(u, √(2π * L / g) / u[2N+2]) # add T to the outcome of the solution
