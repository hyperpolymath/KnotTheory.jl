<!-- SPDX-License-Identifier: MPL-2.0 -->
<!-- Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# Changelog — KnotTheory.jl

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `GaussCode` type for signed Gauss code representation of oriented knot diagrams
- Conversion API: `to_pd`, `to_dt`, `from_dt`, `to_gauss`, `from_gauss`
- Round-trip test coverage for PD↔DT, PD↔Gauss, and braid-word conversions

## [1.0.1]

### Added
- CRG v2 READINESS.md (grade C)
- Deploy dogfood-gate, CRG tests and benchmarks
- EXPLAINME.adoc, TEST-NEEDS.md, PROOF-NEEDS.md

### Changed
- Migrated SCM files to A2ML format in `.machine_readable/6a2/`
- License migration: AGPL-3.0 → MPL-2.0 (source) / MPL-2.0 (LICENSE file)

## [1.0.0]

### Added
- Planar diagram and DT code representations
- Jones, Alexander, Conway, HOMFLY polynomial invariants
- Reidemeister move simplification (R1, R2, R3)
- Seifert matrix + signature computation
- Built-in knot table + lookup
- Braid word conversion utilities
- 1,027 tests passing (unit + E2E + property-based)
