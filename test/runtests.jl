using Test
using KnotTheory

@testset "KnotTheory" begin
    @testset "Basics" begin
        k = unknot()
        @test k isa Knot
        @test k.name == :unknot
        @test crossing_number(k) == 0
    end

    @testset "PD Code" begin
        pd = pdcode([(1, 2, 3, 4, 1)])
        k = Knot(:sample, pd, nothing)
        @test crossing_number(k) == 1
        @test writhe(k) == 1
    end

    @testset "DT Code" begin
        k = trefoil()
        dt = dtcode(k)
        @test dt.code == [4, 6, 2]
        @test crossing_number(k) == 3
    end

    @testset "JSON" begin
        k = Knot(:sample, pdcode([(1, 2, 3, 4, -1)]), DTCode([4, 6, 2]))
        path = joinpath(@__DIR__, "knot.json")
        write_knot_json(path, k)
        k2 = read_knot_json(path)
        @test k2.name == :sample
        @test crossing_number(k2) == 1
        rm(path, force=true)
    end
end
