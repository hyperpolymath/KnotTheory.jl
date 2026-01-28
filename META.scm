;; SPDX-License-Identifier: PMPL-1.0-or-later
;; META.scm - KnotTheory.jl architectural decisions and design rationale
;; Media Type: application/meta+scheme

(define-module (meta knottheory)
  #:use-module (ice-9 match)
  #:export (meta get-adr))

(define meta
  '((metadata
      (version . "0.1.0")
      (schema-version . "1.0.0")
      (created . "2026-01-28")
      (updated . "2026-01-28")
      (project . "KnotTheory.jl")
      (media-type . "application/meta+scheme"))

    (architecture-decisions
      ((adr-001
         (title . "Planar diagram as primary representation")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Knots can be represented via PD codes, DT codes, braid words, or 3D embeddings. Need a representation that balances computational efficiency with algorithmic flexibility.")
         (decision . "Use planar diagram with crossings (4-tuples of arc labels + sign) as primary internal representation. Support PD code and DT code import/export for compatibility. Store components explicitly for multi-component links.")
         (consequences . "Positive: Direct representation for Reidemeister moves, polynomial computation via skein relations. Negative: Not space-efficient for large tables, no direct 3D geometric information."))
       (adr-002
         (title . "Jones polynomial via Kauffman bracket expansion")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Jones polynomial can be computed via skein relations (recursive) or state sum (Kauffman bracket). Need to choose algorithm for clarity and correctness verification.")
         (decision . "Implement Kauffman bracket state sum with smoothing states. Enumerate all 2^n smoothing states, count loops, apply A-polynomial formula, then normalize by writhe for Jones polynomial.")
         (consequences . "Positive: Direct from definition, easy to verify correctness, pedagogically clear. Negative: Exponential complexity O(2^n), not suitable for high-crossing knots without optimization."))
       (adr-003
         (title . "Alexander polynomial via Seifert matrix")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Alexander polynomial can be computed via Fox calculus (knot group), Seifert matrix (from Seifert surface), or skein relations. Need tractable algorithm for initial implementation.")
         (decision . "Use Seifert surface approach: compute Seifert circles by smoothing crossings, build adjacency-based Seifert matrix, compute determinant of (V - V^T t). Use heuristic for small diagrams as scaffold.")
         (consequences . "Positive: Geometric interpretation clear, connects to surface theory. Negative: Heuristic implementation not rigorous, requires refinement for general diagrams, matrix determinant computation needed."))
       (adr-004
         (title . "Reidemeister moves for simplification")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Planar diagrams for the same knot can have different crossing numbers. Need simplification to reduce to minimal crossing representation.")
         (decision . "Implement Reidemeister move I (remove kinks) as r1_simplify. Detect crossings with repeated arc labels and remove them. Future: add R2 and R3 moves for more aggressive simplification.")
         (consequences . "Positive: Reduces diagram complexity, standard topological equivalence. Negative: R1 alone insufficient for minimal crossing number, need R2/R3 for completeness, move detection heuristics needed."))
       (adr-005
         (title . "JSON serialization for knot data exchange")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Need to share knot data with web tools, databases, other languages. Binary formats (JLD2) not portable. PD/DT codes are standard but not self-describing.")
         (decision . "Use JSON3 for read/write with schema: name, pd (list of 5-tuples), dt (list of integers), components (list of arc lists). Human-readable and tool-compatible.")
         (consequences . "Positive: Interoperability, version control friendly, debuggable. Negative: Verbose for large knot tables, no type safety at serialization boundary, requires schema documentation."))
       (adr-006
         (title . "Graphs.jl integration for network analysis")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Knot diagrams have graph structure (arcs as vertices, crossings as edges). May want to analyze connectivity, planarity, cycles. Need to leverage existing graph algorithms.")
         (decision . "Provide to_graph() function converting PlanarDiagram to Graphs.SimpleGraph. Arcs become vertices, crossings define edges. Users can apply standard graph algorithms (shortest paths, connectivity).")
         (consequences . "Positive: Reuse mature graph algorithms, explore graph-theoretic knot invariants. Negative: Graph representation loses crossing sign information, not bijective conversion, interpretation of graph properties unclear."))
       (adr-007
         (title . "Polynomials.jl integration for symbolic manipulation")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Knot invariants are polynomials (Jones, Alexander). May need to manipulate, evaluate, compare polynomials. Need symbolic representation beyond raw coefficient dicts.")
         (decision . "Provide to_polynomial() helper converting Dict{Int,Int} (exponent->coefficient) to Polynomials.Polynomial. Users can evaluate, factor, compare using Polynomials.jl API.")
         (consequences . "Positive: Leverage mature polynomial library, symbolic operations available. Negative: Jones polynomial in Laurent form (negative exponents) requires conversion handling, Polynomials.jl standard polynomial assumed."))
       (adr-008
         (title . "CairoMakie plotting as optional extension")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Visualization useful for teaching and presentations. CairoMakie large dependency (~500MB with artifacts). Not all users need plotting.")
         (decision . "Define plot_pd() as stub in main module. Implement via package extension KnotTheoryCairoMakieExt when CairoMakie loaded. Users opt-in by adding CairoMakie to environment.")
         (consequences . "Positive: No forced dependency, lightweight library by default, users choose visualization backend. Negative: Plotting not available without extension, documentation must clarify opt-in requirement."))
       (adr-009
         (title . "Small knot table with standard examples")
         (status . accepted)
         (date . "2026-01-28")
         (context . "New users need example knots to test algorithms. Full knot tables (KnotInfo) have 10,000+ knots. Want quick access to common examples without large data files.")
         (decision . "Provide knot_table() with unknot, trefoil, figure-eight as initial set. Store DT codes as compact representation. Users import full tables via read_knot_json() if needed.")
         (consequences . "Positive: Fast startup, no data files to manage, covers common teaching examples. Negative: Limited table, users must maintain separate databases for research use, no automatic KnotInfo integration."))
       (adr-010
         (title . "Pure Julia implementation without native dependencies")
         (status . accepted)
         (date . "2026-01-28")
         (context . "Knot invariant computation not performance-critical for small diagrams. Native code (C, Fortran) complicates cross-platform deployment. Want easy installation and wide compatibility.")
         (decision . "Implement all algorithms in pure Julia. Use LinearAlgebra stdlib for matrix operations, JSON3 for serialization, Graphs/Polynomials for integration. Optimize Julia code before considering native extensions.")
         (consequences . "Positive: Easy deployment, cross-platform by default, no build toolchain. Negative: Slower than optimized C for high-crossing knots, no GPU acceleration, but acceptable for typical use (n < 20 crossings).")))

    (development-practices
      (code-style
        (formatter . "julia-format")
        (line-length . 100)
        (naming . "snake_case for functions, PascalCase for types, UPPER_SNAKE for constants")
        (comments . "Docstrings for all exported functions, inline comments for algorithms"))
      (security
        (data-validation . "Check arc labels in range, validate crossing structure")
        (input-sanitization . "Validate JSON schema on deserialization")
        (threat-model . "Assumes trusted input, focus on mathematical correctness"))
      (testing
        (unit-tests . "All exported functions and core algorithms")
        (property-tests . "Polynomial invariant properties (unknot = 1, orientation invariance)")
        (test-vectors . "Known knot invariants from literature (trefoil, figure-eight)")
        (coverage-target . 80))
      (versioning
        (scheme . "SemVer")
        (compatibility . "Julia 1.9+"))
      (documentation
        (api-docs . "Docstrings in source, extracted to docs/")
        (examples . "README with quick start, tutorial notebook in tutorials/")
        (theory . "Algorithm descriptions in docs/README.md")
        (integration . "JSON schema for data exchange"))
      (branching
        (main-branch . "main")
        (feature-branches . "feat/*, fix/*")
        (release-process . "GitHub releases, Julia package registry")))

    (design-rationale
      (why-planar-diagrams
        "Standard representation in knot theory literature"
        "Direct encoding of over/under crossing information"
        "Supports Reidemeister move algorithms naturally"
        "Compatible with PD code and DT code conventions"
        "Enables skein relation computation for polynomials")
      (why-jones-polynomial
        "Most widely used knot invariant in research and teaching"
        "Distinguishes many knots that classical invariants miss"
        "Kauffman bracket provides geometric interpretation"
        "Foundation for quantum topology and TQFT connections"
        "Computable via state sum for moderate crossing numbers")
      (why-alexander-polynomial
        "Classical invariant with rich history (1920s)"
        "Detects knot genus and unknotting information"
        "Connects to knot group and homology theories"
        "Simpler computation than Jones for some diagrams"
        "Educational value for learning knot invariant theory")
      (why-seifert-circles
        "Visual interpretation as smoothed diagram components"
        "Foundation for Seifert surface construction"
        "Relates to knot genus and signature invariants"
        "Simple algorithm: count connected components in smoothing graph"
        "Useful for braid index estimation")
      (why-reidemeister-moves
        "Fundamental theorem: moves generate knot equivalence"
        "Essential for diagram simplification and minimization"
        "Pedagogically important: show topological equivalence"
        "Foundation for invariant verification (move-invariant = knot invariant)"
        "Enables computational topology beyond invariant calculation")
      (why-json-export
        "Enable web-based knot theory tools and databases"
        "Human-readable for debugging and manual inspection"
        "Git-friendly: version control knot data alongside code"
        "Integration with databases (PostgreSQL JSONB, MongoDB)"
        "Cross-language compatibility (Python, JavaScript, R)")
      (why-graphs-integration
        "Knot diagrams are inherently graph-like structures"
        "Leverage mature Graphs.jl algorithms for analysis"
        "Explore graph-theoretic invariants (chromatic polynomial)"
        "Planarity testing, cycle detection, connectivity analysis"
        "Integration with broader Julia graph ecosystem")
      (why-polynomials-integration
        "Invariants are Laurent polynomials in one or more variables"
        "Need symbolic manipulation for comparison and simplification"
        "Polynomial evaluation for specialized invariant values"
        "Factorization to detect product knots"
        "Standard library provides mature polynomial algebra")
      (why-optional-plotting
        "Visualization useful but not core to computation"
        "CairoMakie large dependency, not needed for CLI/server use"
        "Package extensions enable opt-in graphics without bloat"
        "Users can choose backend (Makie, Plots, UnicodePlots)"
        "Keeps library lightweight for batch computation contexts")
      (why-small-knot-table
        "Cover common teaching examples without large data files"
        "Fast startup and testing without database queries"
        "Users import full tables (KnotInfo) via JSON as needed"
        "Separation of computation library from data repository"
        "Encourage users to build domain-specific knot collections")
      (why-pure-julia
        "Knot invariant computation for small diagrams is fast enough"
        "Cross-platform deployment without build complexity"
        "Easy installation via Julia package manager"
        "No binary compatibility issues across OS versions"
        "Julia performance adequate for typical research use (< 20 crossings)")
      (why-dt-and-pd-codes
        "PD code: standard in research, encodes full planar diagram"
        "DT code: compact for alternating knots, classical encoding"
        "Both required for interoperability with existing tools"
        "DT code space-efficient for knot tables"
        "PD code more general, handles non-alternating and links")
      (why-linking-number
        "Classical invariant for 2-component links"
        "Simplest link invariant, pedagogically important"
        "Foundation for more complex link invariants"
        "Efficient computation from planar diagram"
        "Connects to homology and intersection theory"))))

;; Helper function
(define (get-adr id)
  (let ((adrs (assoc-ref meta 'architecture-decisions)))
    (assoc-ref adrs id)))
