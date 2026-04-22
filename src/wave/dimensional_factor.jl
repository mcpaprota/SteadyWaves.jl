module DimensionalFactor

using ..Velocity: VelocityStruct
using ..Wave: WaveStruct

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
    )
end

function dimensional_factor(k,g,rho;L=0,M=0,T=0)
    R = -M
    G = T/2

    K = G - 3R + L

    return k^K * g^G * rho^R
end




# returns compiler that produce factor to multiply dimensional values into dimentionless 
function dimensional_factor_compiler(d,physics)
    g = physics.g
    rho = physics.rho 

    k = (w_c,u) -> w_c.D(w_c,u)/d
    return WaveStruct(
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # eta
	    (w_c, u) -> velocity_struct_factor( k(w_c,u), g),        # v
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # D
	    (w_c, u) -> speed_factor(           k(w_c,u), g),	        # C
	    (w_c, u) -> bernoulli_factor(       k(w_c,u), g),	        # R
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # H
	    (w_c, u) -> speed_factor(           k(w_c,u), g),	        # U
	    (w_c, u) -> flux_factor(            k(w_c,u), g),	        # Q
	    (w_c, u) -> 1,	                                            # N
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # L
	    (w_c, u) -> period_factor(          k(w_c,u), g),	        # T
	    (w_c, u) -> power_factor(           k(w_c,u), g, rho),	    # F
        (w_c, u) -> pressure_factor(        k(w_c,u), g, rho),
        (w_c, u) -> surface_tension_factor( k(w_c,u), g, rho),
	    (w_c, u) -> 1	                                            # raw
    )
end

function dimensional_factor_compiler(L,T, physics)
    g = physics.g
    rho = physics.rho 

    k = (w_c,u) -> 2pi/L
    if L === nothing
        k = (w_c,u) -> w_c.T(w_c,u)/T
    end

    return WaveStruct(
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # eta
	    (w_c, u) -> velocity_struct_factor( k(w_c,u), g),        # v
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # D
	    (w_c, u) -> speed_factor(           k(w_c,u), g),	        # C
	    (w_c, u) -> bernoulli_factor(       k(w_c,u), g),	        # R
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # H
	    (w_c, u) -> speed_factor(           k(w_c,u), g),	        # U
	    (w_c, u) -> flux_factor(            k(w_c,u), g),	        # Q
	    (w_c, u) -> 1,	                                            # N
	    (w_c, u) -> distance_factor(        k(w_c,u)),	            # L
	    (w_c, u) -> period_factor(          k(w_c,u), g),	        # T
	    (w_c, u) -> power_factor(           k(w_c,u), g, rho),	    # F
        (w_c, u) -> pressure_factor(        k(w_c,u), g, rho),
        (w_c, u) -> surface_tension_factor( k(w_c,u), g, rho),
	    (w_c, u) -> 1	                                            # raw
    )
end


end