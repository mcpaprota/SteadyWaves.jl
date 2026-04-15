module DimensionalFactor

using ..Velocity: VelocityStruct
using ..Wave: WaveStruct

function distance_factor(kd,d) 
    return kd/d
end

function period_factor(kd,d,g) 
    return sqrt(g * kd/d)
end

function speed_factor(kd,d,g)
    return sqrt(kd/d / g)
end

function power_factor(kd,d,g,rho)
    return kd/d * sqrt((kd/d / g)^3) / rho
end

function bernoulli_factor(kd,d,g)
    return kd/d / g
end

function flux_factor(kd,d,g)
    return kd/d * sqrt(kd/d/ g) 
end

function pressure_factor(kd,d,g,rho)
    return kd/d / g / rho
end

function surface_tension_factor(kd,d,g,rho)
    return (kd/d)^2 / g / rho
end

function velocity_struct_factor(w_c,u,g,d)
    kd = w_c.D(w_c,u)

    return VelocityStruct(
            speed_factor(kd,d,g),
            speed_factor(kd,d, g),
            bernoulli_factor(kd, d, g),
    )
end

function dimensional_factor(kd,d,g,rho;L=0,M=0,T=0)
    k = kd/d

    R = -M
    G = T/2

    K = G - 3R + L

    return k^K * g^G * rho^R
end

# returns compiler that produce factor to multiply dimensional values into dimentionless 
function dimensional_factor_compiler(d,physics)
    g = physics.g
    rho = physics.rho 

    return WaveStruct(
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # eta
	    (w_c, u) -> velocity_struct_factor(w_c, u, g, d),           # v
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # D
	    (w_c, u) -> speed_factor(      w_c.D(w_c,u), d, g),	        # C
	    (w_c, u) -> bernoulli_factor(  w_c.D(w_c,u), d, g),	        # R
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # H
	    (w_c, u) -> speed_factor(      w_c.D(w_c,u), d, g),	        # U
	    (w_c, u) -> flux_factor(       w_c.D(w_c,u), d, g),	        # Q
	    (w_c, u) -> 1,	                                            # N
	    (w_c, u) -> distance_factor(   w_c.D(w_c,u), d),	        # L
	    (w_c, u) -> period_factor(     w_c.D(w_c,u), d, g),	        # T
	    (w_c, u) -> power_factor(      w_c.D(w_c,u), d, g, rho),	# F
        (w_c, u) -> pressure_factor(   w_c.D(w_c,u), d, g, rho),
        (w_c, u) -> surface_tension_factor( w_c.D(w_c,u), d, g, rho),
	    (w_c, u) -> 1	                                            # raw
    )
end


end