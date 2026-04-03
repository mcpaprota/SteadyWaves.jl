module Linear

using ..Params
using ..Index: IndexStruct
using ..Surface
using ..Wave: WaveStruct
using ..Physics

using NonlinearSolve

function init(d,P,pc,idx, g=G)
    k = Int(pc) == Int(PC_LENGTH) ? 2π / P : linear_wave_number(d, 2π / P, g) # wave number (rad/s)
    u = zeros(idx.U)

    u[idx.D] = k * d

    return k, u
end

function linear_solution(d, P, pc, idx::IndexStruct, compiler, df_compiler; g=G,eta_type::Params.ElevationType =Params.FOURIER_ELEVATION)
    k, u = init(d, P, pc, idx,g)

    df = WaveStruct(u, df_compiler, compiler)

    w = WaveStruct(u, compiler)

    kd = k * d

    amplitudes = [w.D, 0.5 * w.H]

    if eta_type == Params.DIRECT_ELEVATION

        u[idx.eta] = (m -> Surface.fourier_z(amplitudes, w.eta.point_x(m))).(0:idx.N)

    elseif eta_type == Params.FOURIER_ELEVATION

        u[idx.eta[begin:begin+1]] = amplitudes

    end
   

    u[idx.psi[begin]] = 0.5 * w.H / √tanh(kd) # Bk/g
    u[idx.C] = √tanh(kd) # c√(k/g)
    u[idx.D] = w.D # kη̄
    u[idx.Q] = 0 # q√(k³/g)
    u[idx.R] = tanh(kd) / 2 # rk/g
    u[idx.U] = √tanh(kd) # Ū√(k/g)

    return WaveStruct(u,compiler), df
end

function wave_number_condition(k,d,k_0)
    return k * tanh(k *d) - k_0
end


function linear_wave_number_using_solver(d, k_0, ϵ=10^-12)
    system = (du,u,p) -> du[begin] = wave_number_condition(u[begin],d,k_0)

    problem = NonlinearProblem(system, [k_0])

    solution = solve(problem, RobustMultiNewton())

    return solution.u[begin]

end

function linear_wave_number_using_iteration(d,k_0,ϵ=10^-12,max_iter = 1000)
    k = k_0
    i = 0
    while max(abs(k * tanh(k * d) - k_0)) > ϵ 
        k = k_0 / tanh(k * d)

        i += 1
        if i >= max_iter
            break
        end
    end
    return k
end

"""

    linear_wave_number(d, ω, g=G, ϵ=10^-12)

Calculate linear_wave_number `k` based on depth `d`, angular wave frequency `ω`
and gravitational acceleration `g` for given accuracy `ϵ` according to linear wave theory.
"""

function linear_wave_number(d, ω, g=G, ϵ=10^-12)
    k_0 = ω^2 / g # initial guess

    err_k_0 = abs(wave_number_condition(k_0,d,k_0))

    if err_k_0 < ϵ
        return k_0
    end

    k_i = linear_wave_number_using_iteration(d,k_0,ϵ)

    err_k_i = abs(wave_number_condition(k_i,d,k_0))

    if err_k_i < ϵ
        return k_i
    end

    k_s = linear_wave_number_using_solver(d,k_0,ϵ)

    err_k_s = abs(wave_number_condition(k_s,d,k_0))

    min_err = min(err_k_0,err_k_i,err_k_s)

    if min_err == err_k_s
        return k_s
    elseif min_err == err_k_i
        return k_i
    else
        return k_0
    end
end


end