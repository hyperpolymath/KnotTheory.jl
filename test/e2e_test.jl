# SPDX-License-Identifier: MPL-2.0
# E2E pipeline tests for KnotTheory.jl.
# Tests the full pipeline: PD code construction → invariant computation →
# simplification → JSON serialisation → round-trip recovery.

using Test
using KnotTheory
using LinearAlgebra

@testset "E2E Pipeline Tests" begin

    @testset "Full pipeline: PD code → invariants → simplification" begin
        # Step 1: build a trefoil from its PD code.
        t = trefoil()
        @test t isa Knot
        @test crossing_number(t) == 3

        # Step 2: compute all invariants.
        wr = writhe(t)
        @test wr == 3

        n_sc = seifert_circles(t.pd)
        @test n_sc == 2

        V = seifert_matrix(t.pd)
        @test size(V) == (2, 2)

        sig = signature(t.pd)
        @test sig == -2

        det = determinant(t.pd)
        @test det == 3

        alex = alexander_polynomial(t.pd)
        @test abs(sum(values(alex))) == 1

        # Step 3: attempt simplification (trefoil is already minimal).
        reduced = simplify_pd(t.pd)
        @test length(reduced.crossings) == 3
    end

    @testset "Full pipeline: braid word → knot → invariants" begin
        # Build trefoil from braid word and verify full invariant pipeline.
        k = from_braid_word("s1.s1.s1")
        @test crossing_number(k) == 3
        @test writhe(k) == 3
        @test determinant(k.pd) == 3
        @test signature(k.pd) == -2

        # Recover braid word.
        word = to_braid_word(k)
        @test word == "s1.s1.s1"
    end

    @testset "Full pipeline: JSON round-trip preserves all data" begin
        fe = figure_eight()
        path = joinpath(@__DIR__, "e2e_figure_eight.json")
        write_knot_json(path, fe)
        fe2 = read_knot_json(path)
        rm(path, force=true)

        @test fe2.name == :figure_eight
        @test crossing_number(fe2) == 4
        @test fe2.dt !== nothing
        @test fe2.dt.code == [4, 6, 8, 2]
        # Invariants computed from recovered data match originals.
        @test determinant(fe2.pd) == 5
        @test signature(fe2.pd) == 0
    end

    @testset "Full pipeline: unknot has all trivial invariants" begin
        uk = unknot()
        @test crossing_number(uk) == 0
        @test writhe(uk) == 0
        @test seifert_circles(uk.pd) == 0
        @test size(seifert_matrix(uk.pd)) == (0, 0)
        @test signature(uk.pd) == 0
        @test determinant(uk.pd) == 1
        @test alexander_polynomial(uk.pd) == Dict(0 => 1)
        @test conway_polynomial(uk.pd) == Dict(0 => 1)
        @test jones_polynomial(uk.pd; wr=0) == Dict(0 => 1)
    end

    @testset "Error handling: R1 kink removal pipeline" begin
        # A diagram with a kink (R1-removable crossing) should simplify to 0.
        pd = pdcode([(1, 1, 2, 2, 1)])
        @test length(pd.crossings) == 1
        reduced = r1_simplify(pd)
        @test length(reduced.crossings) == 0
        reduced2 = simplify_pd(pd)
        @test length(reduced2.crossings) == 0
    end

    @testset "Error handling: R2 bigon removal pipeline" begin
        pd = pdcode([(1, 3, 2, 4, 1), (2, 4, 5, 6, -1)])
        reduced = r2_simplify(pd)
        @test length(reduced.crossings) < length(pd.crossings)
    end

    @testset "Round-trip consistency: knot table lookup → invariants" begin
        table = knot_table()
        for name in [:trefoil, :figure_eight, Symbol("5_1")]
            entry = lookup_knot(name)
            @test entry !== nothing
            @test entry.crossings == length(entry.dt)
        end
    end

    @testset "Round-trip consistency: graph construction from PD" begin
        t = trefoil()
        g = to_graph(t.pd)
        # Graph must have at least as many vertices as arc labels.
        @test nv(g) >= 6
        @test ne(g) >= 6
    end

end
