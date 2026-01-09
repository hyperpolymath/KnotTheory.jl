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

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```
