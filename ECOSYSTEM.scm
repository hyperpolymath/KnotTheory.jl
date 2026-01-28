;; SPDX-License-Identifier: PMPL-1.0-or-later
;; ECOSYSTEM.scm - Project relationship mapping
;; Media Type: application/vnd.ecosystem+scm

(ecosystem
  (version "1.0")
  (name "KnotTheory.jl")
  (type "computational-mathematics-library")
  (purpose "Julia toolkit for knot theory computation: planar diagrams, polynomial invariants, knot tables, and diagram visualization for topological research and education")

  (position-in-ecosystem
    (role "mathematical-computation-component")
    (layer "application-library")
    (description "Provides computational tools for knot theory within the hyperpolymath mathematical computing ecosystem, supporting research in topology and combinatorics"))

  (related-projects
    ((name . "HackenbushGames.jl")
     (relationship . "sibling-standard")
     (description . "Game theory framework - both work with combinatorial structures and graph representations")
     (integration . "Knot diagrams can be viewed as game positions, braid groups relate to game strategy trees"))
    ((name . "PolyglotFormalisms.jl")
     (relationship . "potential-consumer")
     (description . "Formal proof systems - could use KnotTheory for verified knot invariant proofs")
     (integration . "Formalize knot invariant properties, machine-checked proofs of polynomial calculations"))
    ((name . "Axiom.jl")
     (relationship . "potential-consumer")
     (description . "ML reasoning system - topological data analysis uses knot-theoretic features")
     (integration . "Knot invariants as topological features for ML models, persistent homology integration"))
    ((name . "Causals.jl")
     (relationship . "distant-relation")
     (description . "Causal analysis - both use graph structures but for different purposes")
     (integration . "Graph conversion utilities may be shared, topological properties of causal DAGs"))
    ((name . "KnotInfo")
     (relationship . "inspiration")
     (description . "Comprehensive knot database - KnotTheory provides computational complement")
     (integration . "Compatible data formats for knot tables, import/export knot properties"))
    ((name . "SnapPy")
     (relationship . "inspiration")
     (description . "Python topology library - similar scope but Python-focused")
     (integration . "Compatible PD code and DT code representations for data exchange")))

  (what-this-is
    "A Julia library for computational knot theory with practical algorithms"
    "Planar diagram representation with crossings, arcs, and components"
    "Polynomial invariants: Jones polynomial (Kauffman bracket) and Alexander polynomial"
    "Knot invariants: crossing number, writhe, linking number, Seifert circles"
    "Import/export formats: PD code, DT code, JSON for tool integration"
    "Reidemeister move simplification for diagram reduction"
    "Knot table with lookup for standard knots (unknot, trefoil, figure-eight)"
    "Graph conversion via Graphs.jl for network analysis of knot structure"
    "Polynomial conversion via Polynomials.jl for symbolic manipulation"
    "Optional CairoMakie plotting extension for diagram visualization"
    "Educational scaffold for learning knot theory concepts"
    "Research toolkit for exploring new knot invariants and algorithms")

  (what-this-is-not
    "Not a comprehensive knot database like KnotInfo - focuses on computation"
    "Not a 3D knot visualization tool - uses planar diagram representation"
    "Not a quantum computing library - classical knot invariant algorithms only"
    "Not focused on knot homology theories (Khovanov, Heegaard Floer) - polynomial invariants only"
    "Not a symbolic algebra system - uses numerical computation with Polynomials.jl helpers"
    "Not optimized for massive knot tables - designed for individual knot analysis"
    "Not a topological quantum field theory (TQFT) implementation - practical invariants only"
    "Not a braid group computation tool - estimates braid index but not full braid operations"
    "Not certified for mathematical proof - computational tool for research and exploration"
    "Not a web service - library for local Julia computation"))
