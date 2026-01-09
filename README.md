# KnotTheory.jl

KnotTheory.jl provides a practical Julia toolkit for knot theory: data
structures for planar diagrams, basic invariants, and import/export helpers.
This is an early scaffold intended to grow into a complete library.

## Quick Start

```julia
using KnotTheory

k = trefoil()
println(crossing_number(k))
```

## Scope (initial)

- Planar diagram model with crossings and components.
- Import/export: PD code, DT code, and JSON.
- Basic invariants: crossing number, writhe, linking number.
- Diagram metrics: Seifert circles and braid index estimate.
- Polynomials: Alexander (Seifert matrix heuristic) and Jones (skein).
- Simplification: basic Reidemeister I reduction.
- Example knots: unknot, trefoil, figure-eight.
- Graphs.jl conversion and Polynomials.jl helpers.
- Optional CairoMakie plotting via package extension.

## API Snapshot

```julia
EdgeOrientation, Crossing, PlanarDiagram, DTCode, Knot, Link
pdcode, dtcode, to_dowker
crossing_number, writhe, linking_number
seifert_circles, braid_index_estimate
alexander_polynomial, jones_polynomial
simplify_pd, r1_simplify
knot_table, lookup_knot
write_knot_json, read_knot_json
to_graph, to_polynomial, plot_pd
```

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Docs & Tutorials

- `docs/README.md` for documentation drafts.
- `tutorials/intro.ipynb` for a minimal notebook scaffold.
