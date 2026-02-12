# SPDX-License-Identifier: PMPL-1.0-or-later
using Test
using Graphs
using Polynomials
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
        @test seifert_circles(pd) >= 0
        @test braid_index_estimate(pd) >= 1
    end

    @testset "DT Code" begin
        k = trefoil()
        dt = dtcode(k)
        @test dt.code == [4, 6, 2]
        @test crossing_number(k) == 3
    end

    @testset "Polynomials" begin
        pd = pdcode([(1, 2, 3, 4, 1)])
        alex = alexander_polynomial(pd)
        @test haskey(alex, 0)
        jones = jones_polynomial(pd; wr=1)
        @test !isempty(jones)
        poly, offset = to_polynomial(alex)
        @test poly isa Polynomials.Polynomial
        @test offset isa Int
    end

    @testset "Simplification" begin
        pd = pdcode([(1, 1, 2, 2, 1)])
        reduced = r1_simplify(pd)
        @test length(reduced.crossings) == 0
    end

    @testset "Knot Table" begin
        table = knot_table()
        @test haskey(table, :trefoil)
        @test lookup_knot(:trefoil).crossings == 3
    end

    @testset "Graphs" begin
        pd = pdcode([(1, 2, 3, 4, 1)])
        g = to_graph(pd)
        @test nv(g) >= 4
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
