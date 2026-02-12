;; SPDX-License-Identifier: PMPL-1.0-or-later
;; ECOSYSTEM.scm for KnotTheory.jl

(ecosystem
  (version "1.0")
  (name "KnotTheory.jl")
  (type "julia-package")
  (purpose "Knot theory invariants and planar diagram manipulation")

  (position-in-ecosystem
    (domain "mathematics-topology")
    (role "analytical-library")
    (maturity "alpha")
    (adoption "research-phase"))

  (related-projects
    ((name . "HackenbushGames.jl")
     (relationship . sibling-project)
     (nature . "Combinatorial game theory")))

  (dependencies
    (runtime
      ("Julia" "1.9+")
      ("Polynomials.jl" "4")
      ("Graphs.jl" "1")
      ("JSON3.jl" "1"))))
