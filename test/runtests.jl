using SteadyWaves
using Test

@testset "SteadyWaves.jl" begin
    u = fourier_approx(1, 0.1, 1; cc = 2)
    @test u[22] == u[26]
end
