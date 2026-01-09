# KnotTheory.jl

KnotTheory.jl aims to provide a practical Julia toolkit for knot theory: data
structures for knot diagrams, basic invariants, and utilities for exploration.
This is an early scaffold intended to grow into a complete library.

## Quick Start

```julia
using KnotTheory

k = unknot()
println(k)
```

## Scope (initial)

- Lightweight `Knot` representation and helpers.
- Basic invariants (planned): crossing number, writhe, Alexander/Jones (TBD).
- Diagram utilities (planned): Reidemeister moves and simplification.

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```
