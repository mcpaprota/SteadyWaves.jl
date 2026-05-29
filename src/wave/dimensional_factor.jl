module DimensionalFactor

using ..Velocity: VelocityStruct
using ..Surface: SurfaceStruct, EtaSupportStruct
using ..Wave: WaveStruct
using ..StructOperator: map, combine
using ..Params


function distance_factor(k) 
    return k
end

function period_factor(k,g) 
    return sqrt(g * k)
end

function speed_factor(k,g)
    return sqrt(k / g)
end

function power_factor(k,g,rho)
    return k * sqrt((k / g)^3) / rho
end

function bernoulli_factor(k,g)
    return k / g
end

function flux_factor(k,g)
    return k * sqrt(k/ g) 
end

function pressure_factor(k,g,rho)
    return k / g / rho
end

function surface_tension_factor(k,g,rho)
    return (k)^2 / g / rho
end

function velocity_struct_factor(k,g)
    return VelocityStruct(
            speed_factor(k,g),
            speed_factor(k, g),
            bernoulli_factor(k, g),
            1
    )
end

function surface_struct_factor(k,g,rho)
    func = EtaSupportStruct(
        k,
        k,
        1,
        1/k
    )
    return SurfaceStruct(
        func,
        k,
        k,
        k,
        k,
        (k^4)/(g*rho),
        k,
        func
    )
end


function dimensional_factor(k,g,rho;L=0,M=0,T=0)
    k = k

    R = -M
    G = T/2

    K = G - 3R + L

    return k^K * g^G * rho^R
end

function l_to_k(L)
    return (w_c,u) -> 2pi/L
end

function kd_to_k(d)
    return (w_c,u) -> w_c.D(w_c,u)/d
end

# returns compiler that produce factor to multiply dimensional values into dimentionless 
function dimensional_factor_compiler(d::Number,physics)
    return dimensional_factor_compiler(kd_to_k(d),physics)
end

function dimensional_factor_compiler(definition::Params.Definition,config::Params.ConfigStruct,physics)
    k = nothing
    if config.k_source == Params.K_DEPTH

        k = kd_to_k(definition.d)

    elseif config.k_source == Params.K_PARAMETER

        if definition.L !== nothing

            k = l_to_k(definition.L)
        end
    end

    @assert k !== nothing "Unknown parameters combination"

    return dimensional_factor_compiler(k,physics)
end

function dimensional_factor_compiler(k::Function,physics)
    g = physics.g
    rho = physics.rho 

    return WaveStruct(
	    (w_c, u) -> surface_struct_factor(k(w_c,u),g,rho),	# eta
	    (w_c, u) -> velocity_struct_factor(k(w_c,u),g),           # v
	    (w_c, u) -> distance_factor(   k(w_c,u)),	        # D
	    (w_c, u) -> speed_factor(      k(w_c,u), g),	        # C
        (w_c, u) -> speed_factor(      k(w_c,u), g),	        # c_e
	    (w_c, u) -> bernoulli_factor(  k(w_c,u), g),	        # R
	    (w_c, u) -> distance_factor(   k(w_c,u)),	        # H
	    (w_c, u) -> speed_factor(      k(w_c,u), g),	        # U
	    (w_c, u) -> flux_factor(       k(w_c,u), g),	        # Q
	    (w_c, u) -> 1,	                                            # N
	    (w_c, u) -> distance_factor(   k(w_c,u)),	        # L
	    (w_c, u) -> period_factor(     k(w_c,u), g),	        # T
	    (w_c, u) -> power_factor(      k(w_c,u), g, rho),	# F
        (w_c, u) -> pressure_factor(   k(w_c,u), g, rho),
        (w_c, u) -> surface_tension_factor( k(w_c,u), g, rho),
	    (w_c, u) -> 1	                                            # raw
    )
end

function dimensional(str,df,full_df=nothing)
    if full_df === nothing
        full_df = df
    end

    if typeof(df) <: Number
        if typeof(str) <: Function
            return (args...) -> str((args.*full_df.H)...)/df
        else
            return str/df
        end
    end
    return combine(str,df,(str,df) -> dimensional(str,df,full_df))
end

function dimensionless(str,df,full_df=nothing)
    if full_df === nothing
        full_df = df
    end

    if typeof(df) <: Number
        if typeof(str) <: Function
            return (args...) -> str((args./full_df.H)...)*df
        else
            return str*df
        end
    end
    return combine(str,df,(str,df) -> dimensionless(str,df,full_df))
end

end