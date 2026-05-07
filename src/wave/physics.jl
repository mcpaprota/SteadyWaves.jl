module Physics

export G, RHO, SIGMA

const G = 9.81

const RHO = 1000.0

const SIGMA = 0.0073

struct PhysicsStruct
    g::AbstractFloat
    rho::AbstractFloat
    sigma::AbstractFloat
end

const DEFAULT_PHYSICS = PhysicsStruct(G,RHO,SIGMA)

function validate_parameters(H,L,T,d)

    @assert H > 0

    @assert d > 0

    @assert L === nothing || L > 0 

    @assert T === nothing || T > 0

    @assert (L === nothing) != (T === nothing) 
    
end

function validate_constants(g,rho,sigma)
    
    @assert g > 0

    @assert rho > 0

    @assert sigma > 0
end

function validate_constants(pstruct::PhysicsStruct)

    @assert pstruct.g > 0

    @assert pstruct.rho > 0

    @assert pstruct.sigma > 0
end

end