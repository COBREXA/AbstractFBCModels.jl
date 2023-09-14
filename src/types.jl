
"""
    $(TYPEDEF)

A helper supertype of everything usable as a flux-based constrained linear
metabolic model.

To be functional in most packages that use AbstractFBCModel as an interface,
the following accessors must be overloaded for the concrete model type:
- [`reactions`](@ref)
- [`metabolites`](@ref)
- [`stoichiometry`](@ref)
- [`bounds`](@ref)
- [`objective`](@ref)
"""
abstract type AbstractMetabolicModel end


"""
    $(TYPEDEF)

A nice name for a "nullable" type.
"""
const Maybe{T} = Union{Nothing,T}

"""
    GeneAssociation = Vector{Vector{String}}

An association to genes, represented as a logical formula in a positive
disjunctive normal form (DNF). (The 2nd-level vectors of strings are connected
by "and" to form conjunctions, and the 1st-level vectors of these conjunctions
are connected by "or" to form the DNF.)
"""
const GeneAssociation = Vector{Vector{String}}

"""
    MetaboliteFormula = Dict{String,Int}

Dictionary of atoms and their abundances in a molecule.
"""
const MetaboliteFormula = Dict{String,Int}

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
