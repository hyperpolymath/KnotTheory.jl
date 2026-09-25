# SPDX-License-Identifier: MPL-2.0
# Regression coverage for #42. The normalization fixes in #48 and #52 also
# stabilize the current Alexander-derived HOMFLY output; do not normalize
# away a sign or exponent shift in these comparisons.

@testset "Polynomial crossing-order invariance (#42)" begin
    # Compare both dictionary values and sorted serialization, as downstream
    # descriptor keys do. Dict iteration order itself is not an API guarantee.
    serialise(p) = join([string(k, ":", p[k]) for k in sort!(collect(keys(p)))], ";")
    invariants = (alexander_polynomial, conway_polynomial, homfly_polynomial)

    @testset "Known knots and links" begin
        diagrams = [unknot().pd, trefoil().pd, figure_eight().pd, cinquefoil().pd,
                    from_braid_word("s1.s1").pd,
                    from_braid_word("s1.s1.s2.s2").pd]
        for pd in diagrams
            reordered = PlanarDiagram(reverse(pd.crossings), pd.components)
            for invariant in invariants
                expected = invariant(pd)
                actual = invariant(reordered)
                @test actual == expected
                @test serialise(actual) == serialise(expected)
            end
        end
        # Preserve Conway's normalization, not a positive-leading-coefficient
        # or lowest-degree-zero convention for every polynomial.
        @test alexander_polynomial(figure_eight().pd) == Dict(-1 => -1, 0 => 3, 1 => -1)
        @test conway_polynomial(figure_eight().pd) == Dict(0 => 1, 2 => -1)
        @test alexander_polynomial(from_braid_word("s1.s1").pd) == Dict(0 => -1, 1 => 1)
    end

    # Match quandledb/server/test_br5_fuzz.jl's crossing-order test exactly:
    # seed + 1, strand count, crossing count, then generator BEFORE sign.
    # The older local 60-trial test uses a different seed and draw order.
    rng = MersenneTwister(0xBEAD42 + 1)
    @testset "QuandleDB BR-5 corpus, trial $trial" for trial in 1:200
        nstrands = rand(rng, 3:4)
        ncrossings = rand(rng, 4:12)
        parts = String[]
        for _ in 1:ncrossings
            gen = rand(rng, 1:(nstrands - 1))
            positive = rand(rng, Bool)
            push!(parts, (positive ? "s" : "S") * string(gen))
        end
        pd = from_braid_word(join(parts, ".")).pd
        reordered = PlanarDiagram(shuffle(rng, copy(pd.crossings)), pd.components)
        @testset "$invariant" for invariant in invariants
            expected = invariant(pd)
            actual = invariant(reordered)
            @test actual == expected
            @test serialise(actual) == serialise(expected)
        end

        alex = alexander_polynomial(pd)
        # Pin the representative itself as well as permutation invariance.
        lo, hi = extrema(keys(alex))
        @test lo + hi == (isodd(hi - lo) ? 1 : 0)
        at_one = sum(values(alex))
        @test at_one == 0 ? alex[hi] > 0 : at_one > 0
    end
end
