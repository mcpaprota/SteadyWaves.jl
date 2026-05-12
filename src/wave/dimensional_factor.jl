module DimensionalFactor

using ..Velocity: VelocityStruct
using ..Surface: SurfaceStruct, EtaSupportStruct
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
	    (w_c, u) -> surface_struct_factor(w_c.D(w_c,u)/d,g,rho),	# eta
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

function wrap_eta_support(eta,df_in,df_eta)
     EtaSupportStruct(
       x -> eta.z(x*df_in)/df_eta.z,
       x -> eta.x(x*df_in)/df_eta.x,
       x -> eta.dz_dx_1(x*df_in)/df_eta.dz_dx_1,
       x -> eta.dz_dx_2(x*df_in)/df_eta.dz_dx_2,
     )
end

function dimensional_wave_struct(w,df)
    eta = w.eta

    return WaveStruct(
        SurfaceStruct(
            wrap_eta_support(eta.point, 1, df.eta.point),
            eta.keypoints ./ df.eta.z.z,
            eta.min/df.eta.z.z,
            eta.max/df.eta.z.z,
            eta.avg/df.eta.z.z,
            eta.e_p/df.eta.e_p,
            eta.a ./df.eta.z.z,
            wrap_eta_support(eta.z, df.eta.z.x, df.eta.z)
        ),
        VelocityStruct(
            (x,z) -> w.v.x(x*df.L,z*df.D)/df.v.x,
            (x,z) -> w.v.x(x*df.L,z*df.D)/df.v.x,
            (x,z) -> w.v.psi(x*df.L,z*df.D)/df.v.psi,
            w.v.b,
        ),
        w.D/df.D,
        w.C/df.C,
        w.R/df.R,
        w.H/df.H,
        w.U/df.U,
        w.Q/df.Q,
        w.N,
        w.L/df.L,
        w.T/df.T,
        w.F/df.F,
        (x,z) -> w.P(x*df.L,z*df.D)/df.P,
        w.sigma/df.sigma,
        w.raw
    )

end
end