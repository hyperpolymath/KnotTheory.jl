using Test
using KnotTheory

@testset "KnotTheory" begin
    k = unknot()
    @test k isa Knot
    @test k.name == :unknot
end
