module Physics

export G, RHO

const G = 9.81

const RHO = 1000.0

struct PhysicsStruct
    g::AbstractFloat
    rho::AbstractFloat
end

const DEFAULT_PHYSICS = PhysicsStruct(G,RHO)

end