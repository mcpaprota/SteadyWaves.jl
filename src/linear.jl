module Linear

using ..Params
using ..Index: IndexStruct
using ..Surface
using ..Wave: WaveStruct
using ..Physics


function init(d,P,pc,idx, g=G)
    k = Int(pc) == Int(PC_LENGTH) ? 2π / P : linear_wave_number(d, 2π / P, g) # wave number (rad/s)
    u = zeros(idx.U)

    u[idx.D] = k * d

    return k, u
end

function linear_solution(d, P, pc, idx::IndexStruct, compiler, df_compiler,g=G)
    k, u = init(d, P, pc, idx,g)

    df = WaveStruct(u, df_compiler, compiler)

    w = WaveStruct(u, compiler)

    kd = k * d

    inner = x -> Surface.fourier_z([0.5 * w.H], x, w.D)

    u[idx.eta] = map(inner, (0:idx.N) * π / idx.N) # kη

    u[idx.psi[begin]] = 0.5 * w.H / √tanh(kd) # Bk/g
    u[idx.C] = √tanh(kd) # c√(k/g)
    u[idx.D] = w.D # kη̄
    u[idx.Q] = 0 # q√(k³/g)
    u[idx.R] = tanh(kd) / 2 # rk/g
    u[idx.U] = √tanh(kd) # Ū√(k/g)

    return WaveStruct(u,compiler), df
end

"""

    linear_wave_number(d, ω, g=G, ϵ=10^-12)

Calculate linear_wave_number `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""

function linear_wave_number(d, ω, g=G, ϵ=10^-12)
    k = k₀ = ω^2 / g # initial guess
    while max(abs(k * tanh(k * d) - k₀)) > ϵ
        k = k₀ / tanh(k * d)
    end
    return k
end


end