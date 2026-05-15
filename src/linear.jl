module Linear

using ..Params
using ..Index: Index,IndexStruct
using ..Surface
using ..Wave: Wave, WaveStruct
using ..Physics
using ..DimensionalFactor
using ..Indirect

using NonlinearSolve

function init(d,P,pc,idx, g=G)
    k = Int(pc) == Int(PC_LENGTH) ? 2π / P : linear_wave_number(d, 2π / P, g) # wave number (rad/s)
    u = zeros(idx.U)

    u[idx.D] = k * d

    return k, u
end

function linear_eta_solution(u,w,idx,config)
    amplitudes = [w.D, 0.5 * w.H]

    if config.eta_type == Params.DIRECT_ELEVATION

        u[idx.eta] = (m -> Surface.fourier_z(amplitudes, w.eta.point.x(m))).(0:idx.N)

    elseif config.eta_type == Params.FOURIER_ELEVATION

        u[idx.eta[begin:begin+1]] = amplitudes

    end
end

function linear_angular_frequency(w,config)
    value = config.wave_type == Params.CAPILLARY_WAVE ? 0 : 1
    value += config.wave_type == Params.GRAVITY_WAVE ? 0 : w.sigma

    return value*tanh(w.D)
end
"""
    linear_solution(d,H,P;pc=Params.PC_LENGTH, eta_type::Params.ElevationType=Params.FOURIER_ELEVATION,wave_type::Params.WaveType=Params.GRAVITY_WAVE, g=G, rho=RHO,sigma=SIGMA,N=10)

Approximate solution `u` of a steady wave of height `H` and length `L`
propagating in water of depth `d` using Linear Approximation Method.

# Arguments
- `d`: water depth (m)
- `H`: wave height (m)
- `P`: wave parameter - length `L` (m) or period `T` (s)
- `pc`: parameter criterion; `pc=1`, `pc=PC_LENGTH` - length (default), `pc=2`, `pc=PC_PERIOD` - period
- `eta_type`: elevation representation type;
- `wave_type`: type of forces shaping the wave;
- `N`: number of solution eigenvalues, defaults to `N=10`
- `g`: gravity acceleration (m/s^2), defaults to `g=9.81`
- `rho`: density (kg/m^3), defaults to `ρ=1000`
- `sigma`: surface tension coefficient (kg/s^2), defaults to `σ=0.073`
# Output
- `w`: dimensionless wave struct
- `df`: dimensional factor

"""
function linear_solution(d,H,P; pc=Params.PC_LENGTH, 
    eta_type::Params.ElevationType=Params.FOURIER_ELEVATION,
    wave_type::Params.WaveType=Params.GRAVITY_WAVE,
     g=G, rho=RHO,sigma=SIGMA,N=10)

    physics = Physics.PhysicsStruct(g,rho,sigma)

    Physics.validate_constants(physics)

    idx::IndexStruct = Index.default_indexes(N)

    config = Params.ConfigStruct(
        eta_type = eta_type,
        pc = pc,
    )

    df_compiler = DimensionalFactor.dimensional_factor_compiler(d,physics)

    compiler = WaveStruct(idx,config)


    compiler = Wave.set_compilator_values(
        compiler,
        WaveStruct(H=H,L=Params.L(P,pc),T=Params.T(P,pc)),
        df_compiler
    )

    w, df = linear_solution(d,P,config,idx,compiler,df_compiler,g=physics.g)

    w =  Wave.set_values(w,
            WaveStruct(
            eta = Surface.struct_with_derived_values(w.eta,idx,config.eta_type),
            P = (kx,kz) -> Indirect.indirect_pressure(w,kx,kz),
            F = Indirect.indirect_wave_power(w),
            L = w.L === nothing ? Indirect.indirect_wavelength(w) : nothing,
            T = w.T === nothing ? Indirect.indirect_wave_period(w) : nothing,
        )
    )

    return w, df
end

function linear_solution(d, P, config::Params.ConfigStruct, idx::IndexStruct, compiler, df_compiler; g=G)

    k, u = init(d, P, config.pc, idx,g)

    df = WaveStruct(u, df_compiler, compiler)

    w = WaveStruct(u, compiler)

    linear_eta_solution(u,w,idx,config)

    freq = linear_angular_frequency(w,config)

    omega = √freq #dispersion relation

    u[idx.psi[begin]] = 0.5 * w.H / omega # Bk/g
    u[idx.C] = omega # c√(k/g)
    u[idx.D] = w.D # kη̄
    u[idx.Q] = 0 # q√(k³/g)
    u[idx.R] = freq / 2 # rk/g
    u[idx.U] = omega # Ū√(k/g)

    return WaveStruct(u,compiler), df
end

function dimensionless_linear_solution(config::Params.ConfigStruct, idx::IndexStruct, compiler)
    u = zeros(Index.max_index(idx))

    w = WaveStruct(u, compiler)

    linear_eta_solution(u,w,idx,config)

    freq = linear_angular_frequency(w,config)

    omega = √freq #dispersion relation

    u[idx.psi[begin]] = 0.5 * w.H / omega # Bk/g
    idx.C > 0 ? u[idx.C] = omega :# c√(k/g)
    u[idx.Q] = 0 # q√(k³/g)
    u[idx.R] = freq / 2 # rk/g
    u[idx.U] = omega # Ū√(k/g)

    return WaveStruct(u,compiler), nothing
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

function test_solution_for_linearity(w)
    return sum((w.eta.a[begin+2:end]).^2)/(w.eta.a[begin+1])^2
end

end