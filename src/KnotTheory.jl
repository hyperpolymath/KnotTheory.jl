module KnotTheory

export Knot, unknot

"""
Minimal knot representation.
"""
struct Knot
    name::Symbol
end

"""
Return the unknot.
"""
unknot() = Knot(:unknot)

end # module
