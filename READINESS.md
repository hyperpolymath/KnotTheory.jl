<!--
SPDX-License-Identifier: MPL-2.0
Copyright (c) Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
-->
<!-- Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# Component Readiness — KnotTheory.jl

**Current Grade:** C
**Assessed:** 2026-04-05
**Standard:** [CRG v2.0 STRICT](../standards/component-readiness-grades/)

## Grade rationale (evidence for C)

Works reliably on own project + annotated.

### Evidence

- **Tests:** 1,027 passing (300 unit + 41 E2E + 686 property-based)
- **Annotation:** 122 docstrings across `src/`, EXPLAINME.adoc, TEST-NEEDS.md, PROOF-NEEDS.md, 5 READMEs
- **RSR compliance:** 0-AI-MANIFEST.a2ml, `.machine_readable/6a2/` state files, 14+ workflows, SECURITY/CONTRIBUTING/CODE_OF_CONDUCT
- **Dogfooding:** Consumed by Skein.jl and KRLAdapter.jl as a dependency; 0 regressions across sibling test suites
- **CI:** Clean; panic-attack assail 0 findings

## Gaps preventing higher grades

### Blocks B (6+ diverse external targets)
- No external Julia package beyond the hyperpolymath ecosystem has been tested against this library.
- JuliaHub registration pending.
- No external bug reports have been solicited or received.

### Blocks A (external feedback confirms value)
- Requires B first.

## What to do for B

1. Register on JuliaHub General registry.
2. Solicit 6+ diverse external Julia users working on knot theory, topology,
   or algebraic structures to try the library and report back.
3. Document those trials, issues filed, and fixes shipped.
4. Track the 6 targets explicitly in this file.

## Review cycle

Reassess per release. Next review: on first minor version bump or on any
evidence that test/annotation/RSR claims are no longer accurate.
