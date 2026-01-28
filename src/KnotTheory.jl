# SPDX-License-Identifier: PMPL-1.0-or-later
module KnotTheory

using JSON3
using LinearAlgebra
using Graphs
using Polynomials

export EdgeOrientation, Crossing, PlanarDiagram, DTCode
export Knot, Link
export unknot, trefoil, figure_eight
export crossing_number, writhe, linking_number
export pdcode, dtcode, to_dowker
export write_knot_json, read_knot_json
export seifert_circles, braid_index_estimate
export alexander_polynomial, jones_polynomial
export simplify_pd, r1_simplify
export knot_table, lookup_knot
export to_graph, to_polynomial, plot_pd

@enum EdgeOrientation Over Under

"""
Represents a single crossing with four arc labels and a sign (+1/-1).
"""
struct Crossing
    arcs::NTuple{4, Int}
    sign::Int
end

"""
Planar diagram representation using crossings and optional components.
"""
struct PlanarDiagram
    crossings::Vector{Crossing}
    components::Vector{Vector{Int}}
end

"""
Dowker-Thistlethwaite code (signed even integers).
"""
struct DTCode
    code::Vector{Int}
end

"""
Knot model with optional PD and DT representations.
"""
struct Knot
    name::Symbol
    pd::Union{PlanarDiagram, Nothing}
    dt::Union{DTCode, Nothing}
end

"""
Link model with a planar diagram and component arcs.
"""
struct Link
    name::Symbol
    pd::PlanarDiagram
end

"""
Construct a planar diagram from raw PD tuples.
Each entry is (a, b, c, d, sign).
"""
function pdcode(entries::Vector{NTuple{5, Int}}; components::Vector{Vector{Int}}=Vector{Vector{Int}}())
    crossings = Crossing[]
    for e in entries
        push!(crossings, Crossing((e[1], e[2], e[3], e[4]), e[5]))
    end
    PlanarDiagram(crossings, components)
end

"""
Return the PD code entries (a,b,c,d,sign) for a knot.
"""
function pdcode(knot::Knot)
    knot.pd === nothing && error("knot has no planar diagram")
    [ (c.arcs[1], c.arcs[2], c.arcs[3], c.arcs[4], c.sign) for c in knot.pd.crossings ]
end

"""
Return the DT code for a knot if available or derivable.
"""
function dtcode(knot::Knot)
    if knot.dt !== nothing
        return knot.dt
    end
    knot.pd === nothing && error("knot has no DT code or planar diagram")
    DTCode(to_dowker(knot.pd))
end

"""
Compute Dowker-Thistlethwaite code from a planar diagram.
Assumes a single component and numbered arcs.
"""
function to_dowker(pd::PlanarDiagram)
    mapping = Dict{Int, Int}()
    signmap = Dict{Int, Int}()
    for (idx, crossing) in enumerate(pd.crossings)
        for a in crossing.arcs
            mapping[a] = idx
            signmap[idx] = crossing.sign
        end
    end

    max_arc = isempty(mapping) ? 0 : maximum(keys(mapping))
    code = Int[]
    for odd in 1:2:max_arc
        crossing = mapping[odd]
        even = 0
        for a in pd.crossings[crossing].arcs
            if iseven(a)
                even = a
                break
            end
        end
        if even == 0
            error("could not derive even arc for odd arc $odd")
        end
        push!(code, signmap[crossing] * even)
    end
    code
end

"""
Number of crossings in a knot or link.
"""
function crossing_number(knot::Knot)
    if knot.pd !== nothing
        return length(knot.pd.crossings)
    elseif knot.dt !== nothing
        return length(knot.dt.code)
    end
    0
end

"""
Sum of crossing signs for a knot with a planar diagram.
"""
function writhe(knot::Knot)
    knot.pd === nothing && error("writhe requires a planar diagram")
    sum(c.sign for c in knot.pd.crossings)
end

"""
Compute linking number between two components of a link.
"""
function linking_number(link::Link, comp_a::Int, comp_b::Int)
    comps = link.pd.components
    (comp_a > length(comps) || comp_b > length(comps)) && error("component index out of range")

    comp_map = Dict{Int, Int}()
    for (i, comp) in enumerate(comps)
        for arc in comp
            comp_map[arc] = i
        end
    end

    total = 0
    for c in link.pd.crossings
        arcs = c.arcs
        seen_a = any(get(comp_map, a, 0) == comp_a for a in arcs)
        seen_b = any(get(comp_map, a, 0) == comp_b for a in arcs)
        if seen_a && seen_b
            total += c.sign
        end
    end
    total // 2
end

"""
Count Seifert circles by smoothing crossings.
Uses a simple PD convention and is best for small, well-formed diagrams.
"""
function seifert_circles(pd::PlanarDiagram)
    pairs = Tuple{Int, Int}[]
    for (idx, crossing) in enumerate(pd.crossings)
        a, b, c, d = crossing.arcs
        if crossing.sign >= 0
            push!(pairs, (a, b))
            push!(pairs, (c, d))
        else
            push!(pairs, (b, c))
            push!(pairs, (d, a))
        end
    end

    arc_positions = Dict{Int, Vector{Int}}()
    for (pos, crossing) in enumerate(pd.crossings)
        for arc in crossing.arcs
            push!(get!(arc_positions, arc, Int[]), pos)
        end
    end

    nodes = Int[]
    for (idx, crossing) in enumerate(pd.crossings)
        for _ in 1:4
            push!(nodes, idx)
        end
    end

    adjacency = Dict{Int, Vector{Int}}()
    for (x, y) in pairs
        push!(get!(adjacency, x, Int[]), y)
        push!(get!(adjacency, y, Int[]), x)
    end

    # Count loops in the pairing graph of arcs.
    seen = Set{Int}()
    count = 0
    for arc in keys(adjacency)
        if arc in seen
            continue
        end
        stack = [arc]
        while !isempty(stack)
            cur = pop!(stack)
            if cur in seen
                continue
            end
            push!(seen, cur)
            for n in get(adjacency, cur, Int[])
                if !(n in seen)
                    push!(stack, n)
                end
            end
        end
        count += 1
    end
    count
end

"""
Estimate braid index using Seifert circle count.
"""
braid_index_estimate(pd::PlanarDiagram) = max(1, seifert_circles(pd))

"""
Apply a basic Reidemeister I simplification pass.
Removes crossings with repeated arc labels.
"""
function r1_simplify(pd::PlanarDiagram)
    crossings = Crossing[]
    for c in pd.crossings
        if length(unique(c.arcs)) == 4
            push!(crossings, c)
        end
    end
    PlanarDiagram(crossings, pd.components)
end

"""
Simplify a planar diagram using basic R1 reductions.
"""
simplify_pd(pd::PlanarDiagram) = r1_simplify(pd)

"""
Compute Alexander polynomial from a PD using a simple Seifert matrix heuristic.
Returns a Dict exponent->coefficient in t.
"""
function alexander_polynomial(pd::PlanarDiagram)
    n = seifert_circles(pd)
    if n <= 1
        return Dict(0 => 1)
    end

    # Simple adjacency-based Seifert matrix
    V = zeros(Int, n, n)
    for (i, crossing) in enumerate(pd.crossings)
        a, b, c, d = crossing.arcs

        # Validate arc values are non-negative
        if any(x -> x < 0, [a, b, c, d])
            throw(ArgumentError("Crossing arcs must be non-negative, got $(crossing.arcs)"))
        end

        i1 = (a % n) + 1
        i2 = (c % n) + 1

        # Bounds check
        if i1 < 1 || i1 > n || i2 < 1 || i2 > n
            throw(BoundsError(V, (i1, i2)))
        end

        V[i1, i2] += crossing.sign
    end

    # Compute det(V - V' * t) as polynomial in t
    # Only small matrices are intended.
    t = 1
    size = n
    poly = Dict{Int, Int}()
    for k in 0:size
        poly[k] = 0
    end

    # Use a crude expansion by evaluating determinant at t=0 and t=1
    # as a placeholder; intended for small diagrams.
    det0 = round(Int, det(Matrix{Float64}(V)))
    det1 = round(Int, det(Matrix{Float64}(V - transpose(V))))
    poly[0] = det0
    poly[1] = det1 - det0
    poly
end

const MAX_CROSSINGS_FOR_BRACKET = 20

"""
Compute Jones polynomial via Kauffman bracket expansion.
Returns Dict exponent->coeff where exponent is in quarters of t.
"""
function jones_polynomial(pd::PlanarDiagram; wr::Int=0)
    n = length(pd.crossings)
    if n > MAX_CROSSINGS_FOR_BRACKET
        throw(ArgumentError(
            "Jones polynomial via bracket requires ≤$MAX_CROSSINGS_FOR_BRACKET crossings (got $n)"
        ))
    end
    if n == 0
        return Dict(0 => 1)
    end

    # Build arc positions for pairings.
    arc_positions = Dict{Int, Vector{Int}}()
    for (i, c) in enumerate(pd.crossings)
        for (slot, arc) in enumerate(c.arcs)
            push!(get!(arc_positions, arc, Int[]), 4 * (i - 1) + slot)
        end
    end

    arc_pairs = Tuple{Int, Int}[]
    for positions in values(arc_positions)
        if length(positions) == 2
            push!(arc_pairs, (positions[1], positions[2]))
        end
    end

    function count_loops(pairs::Vector{Tuple{Int, Int}})
        adjacency = Dict{Int, Vector{Int}}()
        for (a, b) in pairs
            push!(get!(adjacency, a, Int[]), b)
            push!(get!(adjacency, b, Int[]), a)
        end
        seen = Set{Int}()
        loops = 0
        for node in keys(adjacency)
            if node in seen
                continue
            end
            stack = [node]
            while !isempty(stack)
                cur = pop!(stack)
                if cur in seen
                    continue
                end
                push!(seen, cur)
                for n in get(adjacency, cur, Int[])
                    if !(n in seen)
                        push!(stack, n)
                    end
                end
            end
            loops += 1
        end
        loops
    end

    memo = Dict{Int, Dict{Int, Int}}()

    function bracket(idx::Int, pairs::Vector{Tuple{Int, Int}})
        if idx > n
            loops = count_loops(pairs)
            # (-A^2 - A^-2)^(loops-1) represented in A powers
            poly = Dict(0 => 1)
            for _ in 1:(loops - 1)
                newpoly = Dict{Int, Int}()
                for (e, c) in poly
                    newpoly[e + 2] = get(newpoly, e + 2, 0) + (-1) * c
                    newpoly[e - 2] = get(newpoly, e - 2, 0) + (-1) * c
                end
                poly = newpoly
            end
            return poly
        end

        # Smoothing for crossing idx
        c = pd.crossings[idx]
        slots = (4 * (idx - 1) + 1, 4 * (idx - 1) + 2, 4 * (idx - 1) + 3, 4 * (idx - 1) + 4)
        a_pairs = vcat(pairs, [(slots[1], slots[2]), (slots[3], slots[4])])
        b_pairs = vcat(pairs, [(slots[2], slots[3]), (slots[4], slots[1])])

        poly_a = bracket(idx + 1, a_pairs)
        poly_b = bracket(idx + 1, b_pairs)

        result = Dict{Int, Int}()
        for (e, c) in poly_a
            result[e + 1] = get(result, e + 1, 0) + c
        end
        for (e, c) in poly_b
            result[e - 1] = get(result, e - 1, 0) + c
        end
        result
    end

    bracket_poly = bracket(1, arc_pairs)

    # Apply writhe normalization: V(t) = (-A)^(-3w) <D>
    a_shift = -3 * wr
    sign = isodd(wr) ? -1 : 1
    jones = Dict{Int, Int}()
    for (e, c) in bracket_poly
        # Convert A^e to t^{-e/4}, track exponent in quarters.
        texp = -(e + a_shift)
        jones[texp] = get(jones, texp, 0) + sign * c
    end
    jones
end

"""
Convert a planar diagram to a Graphs.jl simple graph.
"""
function to_graph(pd::PlanarDiagram)
    max_arc = 0
    for c in pd.crossings
        max_arc = max(max_arc, maximum(c.arcs))
    end
    g = SimpleGraph(max_arc)
    for c in pd.crossings
        add_edge!(g, c.arcs[1], c.arcs[2])
        add_edge!(g, c.arcs[2], c.arcs[3])
        add_edge!(g, c.arcs[3], c.arcs[4])
        add_edge!(g, c.arcs[4], c.arcs[1])
    end
    g
end

"""
Convert a dict exponent->coefficient to a Polynomials.Polynomial.
"""
function to_polynomial(dict::Dict{Int, Int})
    if isempty(dict)
        return Polynomial([0])
    end
    max_exp = maximum(collect(keys(dict)))
    coeffs = zeros(Int, max_exp + 1)
    for (exp, coeff) in dict
        coeffs[exp + 1] = coeff
    end
    Polynomial(coeffs)
end

"""
Plot a planar diagram using CairoMakie if available.
"""
function plot_pd(pd::PlanarDiagram)
    error("plot_pd requires CairoMakie; add it to your environment to enable plotting.")
end

"""
Return the unknot.
"""
unknot() = Knot(:unknot, PlanarDiagram(Crossing[], Vector{Vector{Int}}()), nothing)

"""
Return a trefoil knot (DT code only).
"""
trefoil() = Knot(:trefoil, nothing, DTCode([4, 6, 2]))

"""
Return a figure-eight knot (DT code only).
"""
figure_eight() = Knot(:figure_eight, nothing, DTCode([4, 6, 8, 2]))

"""
Write a knot to JSON at the given path.
"""
function write_knot_json(path::AbstractString, knot::Knot)
    obj = Dict{String, Any}()
    obj["name"] = String(knot.name)
    if knot.pd !== nothing
        obj["pd"] = [ [c.arcs[1], c.arcs[2], c.arcs[3], c.arcs[4], c.sign] for c in knot.pd.crossings ]
        obj["components"] = knot.pd.components
    end
    if knot.dt !== nothing
        obj["dt"] = knot.dt.code
    end
    open(path, "w") do io
        JSON3.write(io, obj)
    end
    nothing
end

"""
Read a knot from JSON produced by write_knot_json.
"""
function read_knot_json(path::AbstractString)
    obj = JSON3.read(read(path, String))
    name = haskey(obj, "name") ? Symbol(String(obj["name"])) : :unnamed

    pd = nothing
    if haskey(obj, "pd")
        entries = Vector{NTuple{5, Int}}()
        for e in obj["pd"]
            push!(entries, (Int(e[1]), Int(e[2]), Int(e[3]), Int(e[4]), Int(e[5])))
        end
        components = Vector{Vector{Int}}()
        if haskey(obj, "components")
            for comp in obj["components"]
                push!(components, [Int(x) for x in comp])
            end
        end
        pd = pdcode(entries; components=components)
    end

    dt = haskey(obj, "dt") ? DTCode([Int(x) for x in obj["dt"]]) : nothing
    Knot(name, pd, dt)
end

"""
Small knot table lookup.
"""
function knot_table()
    Dict(
        :unknot => (name=:unknot, dt=Int[], crossings=0),
        :trefoil => (name=:trefoil, dt=[4, 6, 2], crossings=3),
        :figure_eight => (name=:figure_eight, dt=[4, 6, 8, 2], crossings=4),
    )
end

"""
Lookup a knot entry by name.
"""
lookup_knot(name::Symbol) = get(knot_table(), name, nothing)

end # module
