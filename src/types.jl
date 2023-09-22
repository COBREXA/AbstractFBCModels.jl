
using SparseArrays

"""
    $(TYPEDEF)

A helper supertype of flux balanace based metabolic model.

To create concrete subtypes, it is required to implement all the functions
listed in the module definition.
"""
abstract type AbstractFBCModel end

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

"""
A shortname for a sparse matrix.
"""
const SparseMat = SparseMatrixCSC{Float64,Int}

"""
A shortname for a sparse vector.
"""
const SparseVec = SparseVector{Float64,Int}

"""
    GeneAssociationsDNF = Vector{Vector{String}}

Disjunctive normal form of simple gene associations. For example, `[[A, B],
[B]]` represents two isozymes where the first requires both genes `A` and `B`,
while the second isozyme only requires gene `C`.

This string representation is typically used to represent gene reaction rules,
but does not contain any subunit stoichiometry of kinetic information of the
isozymes. See [`Isozyme`}(@ref) for a more complete structure.
"""
const GeneAssociationsDNF = Vector{Vector{String}}

"""
    Annotations = Dict{String,Vector{String}}

Dictionary used to store (possible multiple) standardized annotations of
something, such as a [`Metabolite`](@ref) and a [`Reaction`](@ref).

# Example
```
Annotations("PubChem" => ["CID12345", "CID54321"])
```
"""
const Annotations = Dict{String,Vector{String}}

"""
    Notes = Dict{String,Vector{String}}

Free-form notes about something (e.g. a [`Gene`](@ref)), categorized by
"topic".
"""
const Notes = Dict{String,Vector{String}}
