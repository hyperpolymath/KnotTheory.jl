# SPDX-License-Identifier: MPL-2.0
# (PMPL-1.0-or-later preferred; MPL-2.0 required for Julia ecosystem)
# BenchmarkTools benchmarks for KnotTheory.jl.

using BenchmarkTools
using KnotTheory

println("=== KnotTheory.jl Benchmarks ===")

# --- Knot construction ---

println("\n-- Knot construction --")

b_unknot   = @benchmark unknot()
b_trefoil  = @benchmark trefoil()
b_fig8     = @benchmark figure_eight()
b_cinq     = @benchmark cinquefoil()
println("unknot():      ", median(b_unknot))
println("trefoil():     ", median(b_trefoil))
println("figure_eight(): ", median(b_fig8))
println("cinquefoil():  ", median(b_cinq))

# --- Invariant computation ---

println("\n-- Invariant computation (trefoil) --")

t = trefoil()
fe = figure_eight()

b_writhe     = @benchmark writhe($t)
b_seifert_c  = @benchmark seifert_circles($t.pd)
b_seifert_m  = @benchmark seifert_matrix($t.pd)
b_signature  = @benchmark signature($t.pd)
b_determinant = @benchmark determinant($t.pd)
b_alexander  = @benchmark alexander_polynomial($t.pd)
b_conway     = @benchmark conway_polynomial($t.pd)
b_jones      = @benchmark jones_polynomial($t.pd; wr=writhe($t))
b_homfly     = @benchmark homfly_polynomial($t.pd)
println("writhe:               ", median(b_writhe))
println("seifert_circles:      ", median(b_seifert_c))
println("seifert_matrix:       ", median(b_seifert_m))
println("signature:            ", median(b_signature))
println("determinant:          ", median(b_determinant))
println("alexander_polynomial: ", median(b_alexander))
println("conway_polynomial:    ", median(b_conway))
println("jones_polynomial:     ", median(b_jones))
println("homfly_polynomial:    ", median(b_homfly))

# --- Simplification ---

println("\n-- Simplification --")

kink_pd = pdcode([(1, 1, 2, 2, 1)])
b_r1 = @benchmark r1_simplify($kink_pd)
b_simplify_trefoil = @benchmark simplify_pd($t.pd)
b_simplify_fe = @benchmark simplify_pd($fe.pd)
println("r1_simplify (kink):       ", median(b_r1))
println("simplify_pd (trefoil):    ", median(b_simplify_trefoil))
println("simplify_pd (fig. eight): ", median(b_simplify_fe))

# --- JSON round-trip ---

println("\n-- JSON round-trip --")

import Serialization, Mmap
tmp_path = joinpath(tempdir(), "bench_knot.json")

b_write = @benchmark write_knot_json($tmp_path, $t)
println("write_knot_json (trefoil): ", median(b_write))

b_read = @benchmark read_knot_json($tmp_path)
println("read_knot_json (trefoil):  ", median(b_read))
rm(tmp_path, force=true)
