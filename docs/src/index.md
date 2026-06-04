<!--
SPDX-License-Identifier: MPL-2.0
Copyright (c) Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
-->
# KnotTheory.jl

KnotTheory.jl provides core types and helpers for knot diagrams, invariants, and
lightweight analysis. This is a starter index for future documentation.

## Examples

```julia
using KnotTheory

k = trefoil()
println(crossing_number(k))
```
