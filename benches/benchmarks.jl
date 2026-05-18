# SPDX-License-Identifier: MPL-2.0
# (MPL-2.0 is automatic legal fallback until PMPL is formally recognised)
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# KnotTheory.jl — BenchmarkTools benchmark suite with Six Sigma classification.
#
# Run:
#   julia --project=. benches/benchmarks.jl
#
# Six Sigma tiers (relative to BASELINES below):
#   UNACCEPTABLE  : >50 % regression  — hard fail
#   ACCEPTABLE    : 20–50 % regression — soft fail
#   ORDINARY      : ±20 %             — pass
#   EXTRAORDINARY : >20 % improvement — pass + flag
#
# First run: BASELINES is empty.  The script prints "[BASELINE]" lines.
# Copy the printed ns values into BASELINES for subsequent CI comparisons.

using BenchmarkTools
using KnotTheory

# ── Six Sigma baselines (nanoseconds) ─────────────────────────────────────────
# Populate from a "[BASELINE]" run.  Leave 0.0 to treat a key as unset.
const BASELINES = Dict{String, Float64}(
    # construction
    "unknot()"         => 0.0,
    "trefoil()"        => 0.0,
    "figure_eight()"   => 0.0,
    "cinquefoil()"     => 0.0,
    # invariants (trefoil)
    "writhe"                => 0.0,
    "seifert_circles"       => 0.0,
    "seifert_matrix"        => 0.0,
    "signature"             => 0.0,
    "determinant"           => 0.0,
    "alexander_polynomial"  => 0.0,
    "conway_polynomial"     => 0.0,
    "jones_polynomial"      => 0.0,
    "homfly_polynomial"     => 0.0,
    # simplification
    "r1_simplify (kink)"        => 0.0,
    "simplify_pd (trefoil)"     => 0.0,
    "simplify_pd (fig. eight)"  => 0.0,
    # JSON round-trip
    "write_knot_json"  => 0.0,
    "read_knot_json"   => 0.0,
)

# ── Six Sigma classifier ───────────────────────────────────────────────────────
const _SIGMA_COUNTS = Dict(:baseline => 0, :extraordinary => 0,
                            :ordinary => 0, :acceptable => 0, :unacceptable => 0)

function classify_sigma(label::String, measured_ns::Float64)::Symbol
    baseline = get(BASELINES, label, 0.0)
    if baseline == 0.0
        @printf("  [BASELINE]      %-38s  %.1f ns\n", label, measured_ns)
        _SIGMA_COUNTS[:baseline] += 1
        return :baseline
    end
    pct = (measured_ns - baseline) / baseline * 100.0
    if pct > 50.0
        @printf("  [UNACCEPTABLE]  %-38s  %+.1f %%  HARD FAIL\n", label, pct)
        _SIGMA_COUNTS[:unacceptable] += 1
        return :unacceptable
    elseif pct > 20.0
        @printf("  [ACCEPTABLE]    %-38s  %+.1f %%  soft fail\n", label, pct)
        _SIGMA_COUNTS[:acceptable] += 1
        return :acceptable
    elseif pct >= -20.0
        @printf("  [ORDINARY]      %-38s  %+.1f %%\n", label, pct)
        _SIGMA_COUNTS[:ordinary] += 1
        return :ordinary
    else
        @printf("  [EXTRAORDINARY] %-38s  %+.1f %%  improvement\n", label, pct)
        _SIGMA_COUNTS[:extraordinary] += 1
        return :extraordinary
    end
end

function bench_ns(b::BenchmarkTools.Trial)::Float64
    time(median(b))
end

# ── Helpers ────────────────────────────────────────────────────────────────────
using Printf

println("=== KnotTheory.jl Benchmarks (Six Sigma) ===")
println("Date: $(Dates.today())\n")
using Dates

# ── 1. Knot construction ───────────────────────────────────────────────────────
println("─── Knot construction ─────────────────────────────────────────────────")

b_unknot  = @benchmark unknot()
b_trefoil = @benchmark trefoil()
b_fig8    = @benchmark figure_eight()
b_cinq    = @benchmark cinquefoil()

classify_sigma("unknot()",       bench_ns(b_unknot))
classify_sigma("trefoil()",      bench_ns(b_trefoil))
classify_sigma("figure_eight()", bench_ns(b_fig8))
classify_sigma("cinquefoil()",   bench_ns(b_cinq))

println()

# ── 2. Invariant computation (trefoil) ─────────────────────────────────────────
println("─── Invariant computation (trefoil) ───────────────────────────────────")

t  = trefoil()
fe = figure_eight()
wr = writhe(t)

b_writhe     = @benchmark writhe($t)
b_seifert_c  = @benchmark seifert_circles($t.pd)
b_seifert_m  = @benchmark seifert_matrix($t.pd)
b_signature  = @benchmark signature($t.pd)
b_determinant = @benchmark determinant($t.pd)
b_alexander  = @benchmark alexander_polynomial($t.pd)
b_conway     = @benchmark conway_polynomial($t.pd)
b_jones      = @benchmark jones_polynomial($t.pd; wr=$wr)
b_homfly     = @benchmark homfly_polynomial($t.pd)

classify_sigma("writhe",               bench_ns(b_writhe))
classify_sigma("seifert_circles",      bench_ns(b_seifert_c))
classify_sigma("seifert_matrix",       bench_ns(b_seifert_m))
classify_sigma("signature",            bench_ns(b_signature))
classify_sigma("determinant",          bench_ns(b_determinant))
classify_sigma("alexander_polynomial", bench_ns(b_alexander))
classify_sigma("conway_polynomial",    bench_ns(b_conway))
classify_sigma("jones_polynomial",     bench_ns(b_jones))
classify_sigma("homfly_polynomial",    bench_ns(b_homfly))

println()

# ── 3. Simplification ──────────────────────────────────────────────────────────
println("─── Simplification ────────────────────────────────────────────────────")

kink_pd = pdcode([(1, 1, 2, 2, 1)])

b_r1            = @benchmark r1_simplify($kink_pd)
b_simplify_t    = @benchmark simplify_pd($t.pd)
b_simplify_fe   = @benchmark simplify_pd($fe.pd)

classify_sigma("r1_simplify (kink)",       bench_ns(b_r1))
classify_sigma("simplify_pd (trefoil)",    bench_ns(b_simplify_t))
classify_sigma("simplify_pd (fig. eight)", bench_ns(b_simplify_fe))

println()

# ── 4. JSON round-trip ─────────────────────────────────────────────────────────
println("─── JSON round-trip ───────────────────────────────────────────────────")

tmp_path = joinpath(tempdir(), "bench_knot_sigma.json")

b_write = @benchmark write_knot_json($tmp_path, $t)
classify_sigma("write_knot_json", bench_ns(b_write))

b_read = @benchmark read_knot_json($tmp_path)
classify_sigma("read_knot_json", bench_ns(b_read))

rm(tmp_path, force=true)

println()

# ── Summary ────────────────────────────────────────────────────────────────────
println("─── Six Sigma Summary ─────────────────────────────────────────────────")
total = sum(values(_SIGMA_COUNTS))
if _SIGMA_COUNTS[:baseline] == total
    println("  BASELINE RUN — no prior measurements.  Record the ns values above.")
    println("  Copy them into the BASELINES dict and re-run to classify.")
else
    hard_fails = _SIGMA_COUNTS[:unacceptable]
    soft_fails = _SIGMA_COUNTS[:acceptable]
    @printf("  Baseline:      %d\n", _SIGMA_COUNTS[:baseline])
    @printf("  Extraordinary: %d\n", _SIGMA_COUNTS[:extraordinary])
    @printf("  Ordinary:      %d\n", _SIGMA_COUNTS[:ordinary])
    @printf("  Acceptable:    %d  (soft fail)\n", soft_fails)
    @printf("  Unacceptable:  %d  (HARD FAIL)\n", hard_fails)
    println()
    if hard_fails > 0
        println("  RESULT: FAIL — $(hard_fails) hard regression(s)")
    elseif soft_fails > 0
        println("  RESULT: WARN — $(soft_fails) soft regression(s), no hard fails")
    else
        println("  RESULT: PASS")
    end
end

println("\n=== Done ===")
