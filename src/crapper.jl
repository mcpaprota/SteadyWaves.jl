module Crapper

using ..Physics: PhysicsStruct
using ..Wave: WaveStruct


function c0_c_ratio_squared(;steepness)
    return  sqrt(1 + (π^2 * steepness^2) / 4)
end 

function c0_c_ratio_squared(;kH)
    return  sqrt(1 + (kH^2) / 16)
end 

function c0_squared(k, physics::PhysicsStruct)
    return k * k / physics.g * physics.sigma/physics.rho
end

function c0_squared(w::WaveStruct)
    return w.sigma
end

function crapper_test_of_solution(w)
    return w.C^2 * c0_c_ratio_squared(kH=w.H) - c0_squared(w)
end

end