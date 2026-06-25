<!--
SPDX-License-Identifier: CC-BY-SA-4.0
SPDX-FileCopyrightText: 2025-2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
-->

[![Project Topology](https://img.shields.io/badge/Project-Topology-9558B2)](TOPOLOGY.md)
[![Completion Status](https://img.shields.io/badge/Completion-100%25-brightgreen)](TOPOLOGY.md) [![OpenSSF Best Practices](https://img.shields.io/badge/OpenSSF-Best_Practices-green?logo=opensourcesecurity)](https://www.bestpractices.dev/en/projects/new?repo_url=https://github.com/hyperpolymath/KnotTheory.jl)
[![License](https://img.shields.io/badge/License-MPL--2.0-blue.svg)](LICENSE) <embed
src="https://api.thegreenwebfoundation.org/greencheckimage/github.com"
data-link="https://www.thegreenwebfoundation.org/green-web-check/?url=github.com" />

A comprehensive Julia toolkit for computational knot theory: planar
diagram data structures, classical invariants, polynomial invariants,
Seifert theory, Reidemeister simplification, braid word interop, and
import/export helpers.

# Installation

## From Julia REPL

```julia
using Pkg
Pkg.add("KnotTheory")
```

## From Git (Development)

```julia
using Pkg
Pkg.add(url="https://github.com/hyperpolymath/KnotTheory.jl")
```

# Quick Start

```julia
using KnotTheory

k = trefoil()
println(crossing_number(k))       # 3
println(alexander_polynomial(k))   # t^-1 - 1 + t
println(jones_polynomial(k))       # -t^-4 + t^-3 + t^-1
```

# Features

- **Planar diagram model** with oriented crossings and multi-component
  links.

- **Code representations**: PD code, DT/Dowker code, signed Gauss code,
  JSON serialization.

- **Classical invariants**: crossing number, writhe, linking number,
  signature, determinant.

- **Polynomial invariants**: Alexander, Jones, Conway, and HOMFLY-PT.

- **Seifert theory**: Seifert circles, Seifert matrix, braid index
  estimate.

- **Reidemeister simplification**: R1, R2, R3 moves and combined
  simplifier.

- **Braid word interop**: convert between planar diagrams and braid
  words (TANGLE compatibility).

- **Knot table**: built-in catalogue with named knots up to standard
  tables.

- **Graph conversion**: Graphs.jl integration via `to_graph`.

- **Polynomial helpers**: Polynomials.jl conversion via `to_polynomial`.

- **Optional plotting**: CairoMakie-based diagram rendering via package
  extension.

# API Reference

## Types

| Type              | Description                                         |
|-------------------|-----------------------------------------------------|
| `EdgeOrientation` | Enum for edge direction (`Over`, `Under`)           |
| `Crossing`        | Single crossing with strand indices and orientation |
| `PlanarDiagram`   | Full planar diagram with crossings and components   |
| `DTCode`          | Dowker-Thistlethwaite code representation           |
| `GaussCode`       | Signed Gauss code representation                    |
| `Knot`            | Named knot wrapper (e.g. `trefoil()`)               |
| `Link`            | Named link wrapper for multi-component objects      |

## Constructors

| Function            | Description                      |
|---------------------|----------------------------------|
| `unknot()`          | The unknot (zero crossings)      |
| `trefoil()`         | Trefoil knot (3_1)               |
| `figure_eight()`    | Figure-eight knot (4_1)          |
| `cinquefoil()`      | Cinquefoil knot (5_1)            |
| `knot_table()`      | Return the built-in knot table   |
| `lookup_knot(name)` | Look up knot table entry by name |

## Classical Invariants

| Function             | Description                            |
|----------------------|----------------------------------------|
| `crossing_number(k)` | Minimum crossing number                |
| `writhe(pd)`         | Sum of crossing signs                  |
| `linking_number(pd)` | Linking number for two-component links |
| `signature(k)`       | Knot signature (from Seifert matrix)   |
| `determinant(k)`     | Knot determinant (                     |
| det(V + V^T)         | )                                      |

## Polynomial Invariants

| Function                  | Description                                     |
|---------------------------|-------------------------------------------------|
| `alexander_polynomial(k)` | Alexander polynomial via Seifert matrix         |
| `jones_polynomial(k)`     | Jones polynomial via skein relation             |
| `conway_polynomial(k)`    | Conway polynomial (substitution from Alexander) |
| `homfly_polynomial(k)`    | HOMFLY-PT two-variable polynomial               |

## Seifert Theory

| Function | Description |
|----|----|
| `seifert_circles(pd)` | Seifert circle decomposition |
| `seifert_circles_with_map(pd)` | Seifert circles with strand-to-circle mapping |
| `seifert_matrix(pd)` | Seifert matrix computation |
| `braid_index_estimate(pd)` | Lower bound on braid index from Seifert circles |

## Conversion & Serialization

| Function | Description |
|----|----|
| `pdcode(k)` | Planar diagram code (list of crossing tuples) |
| `dtcode(k)` | Dowker-Thistlethwaite code |
| `to_pd(x)` | Canonical conversion to planar diagram |
| `to_dt(x)` | Convert to DT code |
| `from_dt(dt)` | Convert DT code to planar diagram |
| `to_gauss(x)` | Convert to signed Gauss code |
| `from_gauss(g)` | Convert signed Gauss code to planar diagram |
| `write_knot_json(file,` `k)` | Serialize knot data to JSON |
| `read_knot_json(file)` | Deserialize knot data from JSON |

## Simplification

| Function          | Description                                |
|-------------------|--------------------------------------------|
| `simplify_pd(pd)` | Apply all Reidemeister moves until stable  |
| `r1_simplify(pd)` | Reidemeister I: remove kinks               |
| `r2_simplify(pd)` | Reidemeister II: cancel opposing crossings |
| `r3_simplify(pd)` | Reidemeister III: triangle move            |

## Braid Words (TANGLE Interop)

| Function                | Description                                  |
|-------------------------|----------------------------------------------|
| `from_braid_word(word)` | Construct a knot from braid word             |
| `to_braid_word(k)`      | Convert knot or planar diagram to braid word |

## Utilities

| Function              | Description                                    |
|-----------------------|------------------------------------------------|
| `to_graph(pd)`        | Convert to Graphs.jl graph structure           |
| `to_polynomial(expr)` | Convert to Polynomials.jl polynomial           |
| `plot_pd(pd)`         | Render diagram (requires CairoMakie extension) |

# Development

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

285 tests across 25 test sets covering all exported functions, invariant
consistency, and known values from knot tables.

# External Integration API

For external consumers (e.g. `Skein.jl`), the minimal surface is:

- `PlanarDiagram`, `Crossing`, `DTCode`, `GaussCode` (core
  representation types)

- `to_pd`, `to_dt`, `to_gauss`, `from_dt`, `from_gauss` (pure
  conversions)

- `crossing_number`, `writhe`, `seifert_circles`, `seifert_matrix`

- `alexander_polynomial`, `jones_polynomial`, `conway_polynomial`,
  `homfly_polynomial`

- `signature`, `determinant`

- `simplify_pd`, `r1_simplify`, `r2_simplify`, `r3_simplify`

This package intentionally does not include persistence, indexing, or
database logic.

# Docs & Tutorials

- `docs/README.md` for documentation drafts.

- `tutorials/intro.ipynb` for a minimal notebook scaffold.

# References & Bibliography

## Textbooks

- Adams, C.C. *The Knot Book*. American Mathematical Society, 2004.

- Lickorish, W.B.R. *An Introduction to Knot Theory*. Springer, 1997.

- Rolfsen, D. *Knots and Links*. AMS Chelsea Publishing, 1976/2003.

- Murasugi, K. *Knot Theory and Its Applications*. Birkhauser, 1996.

- Kauffman, L.H. *Knots and Physics*. 3rd ed., World Scientific, 2001.

- Cromwell, P.R. *Knots and Links*. Cambridge University Press, 2004.

## Key Papers

- Fox, R.H. "Free differential calculus. I." *Annals of Mathematics*
  57(3), 1953, pp. 547–560.

- Freyd, P. et al. "A new polynomial invariant of knots and links."
  *Bulletin of the AMS* 12(2), 1985, pp. 239–246.

- Seifert, H. "Uber das Geschlecht von Knoten." *Mathematische Annalen*
  110, 1935, pp. 571–592.

- Jones, V.F.R. "A polynomial invariant for knots via von Neumann
  algebras." *Bulletin of the AMS* 12(1), 1985, pp. 103–111.

# License

SPDX-License-Identifier: CC-BY-SA-4.0 See [LICENSE](LICENSE).
