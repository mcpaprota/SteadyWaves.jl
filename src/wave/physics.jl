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

end