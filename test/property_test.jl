# SPDX-License-Identifier: MPL-2.0
# Property-based invariant tests for KnotTheory.jl.

using Test
using KnotTheory
using LinearAlgebra

@testset "Property-Based Tests" begin

    @testset "Alexander polynomial is symmetric for all standard knots" begin
        knot_fns = [trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            alex = alexander_polynomial(k.pd)
            for (e, c) in alex
                @test get(alex, -e, 0) == c
            end
        end
    end

    @testset "Alexander polynomial evaluates to ±1 at t=1 for all standard knots" begin
        knot_fns = [trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            alex = alexander_polynomial(k.pd)
            @test abs(sum(values(alex))) == 1
        end
    end

    @testset "Seifert matrix dimension formula: dim = c - s + 1" begin
        knot_fns = [trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            c = length(k.pd.crossings)
            s = seifert_circles(k.pd)
            V = seifert_matrix(k.pd)
            expected_dim = c - s + 1
            @test size(V, 1) == expected_dim
            @test size(V, 2) == expected_dim
        end
    end

    @testset "Determinant equals |det(V + V')| for all standard knots" begin
        knot_fns = [trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            V = seifert_matrix(k.pd)
            d = determinant(k.pd)
            computed = abs(round(Int, det(Float64.(V + transpose(V)))))
            @test d == computed
        end
    end

    @testset "Signature is always even (knot invariant)" begin
        knot_fns = [unknot, trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            sig = signature(k.pd)
            @test iseven(sig)
        end
    end

    @testset "Conway polynomial has only even powers for knots" begin
        knot_fns = [trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            conway = conway_polynomial(k.pd)
            for (e, _) in conway
                @test iseven(e)
            end
        end
    end

    @testset "crossing_number always equals PD crossing count" begin
        knot_fns = [unknot, trefoil, figure_eight, cinquefoil]
        for _ in 1:50
            k = rand(knot_fns)()
            @test crossing_number(k) == length(k.pd.crossings)
        end
    end

    @testset "simplify_pd is idempotent on minimal diagrams" begin
        knot_fns = [trefoil, figure_eight]
        for _ in 1:20
            k = rand(knot_fns)()
            r1 = simplify_pd(k.pd)
            r2 = simplify_pd(r1)
            @test length(r1.crossings) == length(r2.crossings)
        end
    end

    @testset "lookup_knot returns nothing for non-existent names" begin
        bogus_names = [Symbol("x_9999"), :not_a_knot, Symbol("bogus_$(rand(1:1000))")]
        for _ in 1:50
            name = rand(bogus_names)
            @test lookup_knot(name) === nothing
        end
    end

    @testset "knot_table entries: DT code length equals crossing number" begin
        table = knot_table()
        for _ in 1:30
            entry = rand(values(table))
            @test length(entry.dt) == entry.crossings
        end
    end

end
