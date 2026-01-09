module KnotTheory

using JSON3

export EdgeOrientation, Crossing, PlanarDiagram, DTCode
export Knot, Link
export unknot, trefoil, figure_eight
export crossing_number, writhe, linking_number
export pdcode, dtcode, to_dowker
export write_knot_json, read_knot_json

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

end # module
