# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Sealed-container packaging for KnotTheory.jl — the estate package policy's
# sanctioned escape hatch (0-canon/rsr/3-practice/LANGUAGE-POLICY.adoc,
# RULED 2026-05-18: Guix primary + sealed-container escape; NO Nix mirror).
#
# Guix is the estate primary, but this package depends on the unregistered
# estate package AcceleratorGate.jl as a sibling path dependency, which no
# Guix channel carries — exactly the not-in-Guix tail the escape hatch is
# for. Build and self-test with Podman:
#
#   podman build -t knottheory .
#   podman run --rm knottheory
#
# The steps mirror .github/workflows/ci.yml, which is the green reference
# environment for this repo.

FROM docker.io/julia:1.10-bookworm

WORKDIR /estate/KnotTheory.jl
COPY . .

# KnotTheory depends on the unregistered estate package AcceleratorGate.jl as
# a path-dep (Manifest: ../AcceleratorGate.jl) — clone it as a sibling.
RUN git clone --depth 1 https://github.com/hyperpolymath/AcceleratorGate.jl.git /estate/AcceleratorGate.jl

# The committed Manifest.toml is resolved with a specific Julia version, so
# discard it and re-resolve fresh for this image's Julia, devving the
# path-dep (does not alter Project.toml).
RUN rm -f Manifest.toml \
    && julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path="/estate/AcceleratorGate.jl")); Pkg.instantiate()'

# Default: run the full test suite — the same matrix ci.yml runs on Julia
# 1.10 and 1.11.
ENTRYPOINT ["julia", "--project=.", "-e", "using Pkg; Pkg.test()"]
