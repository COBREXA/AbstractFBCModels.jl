
using SparseArrays

"""
$(TYPEDEF)

A supertype of all Flux-balance-Based Constrained metabolic models that share
the common interface defined by `AbstractFBCModels.jl`.

To use the interface for your type, make it a subtype of
[`AbstractFBCModel`](@ref) and provide methods for various functions used with
the model. Use [`accessors`](@ref) to find the current list of methods, and
utilize [`run_fbcmodel_type_tests`](@ref) and [`run_fbcmodel_file_tests`](@ref)
to test the completeness and compatibility of your implementation with the
assumptions of the interface.
"""
abstract type AbstractFBCModel end

"""
$(TYPEDEF)

A nice name for a "nullable" type.
"""
const Maybe{T} = Union{Nothing,T}

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
    GeneAssociationDNF = Vector{Vector{String}}

Disjunctive normal form of simple gene associations. For example, `[[A, B],
[B]]` represents two possibilities to run a given reaction, where the first
requires both gene products `A` and `B`, while the second possibility only
requires gene product `C`.
"""
const GeneAssociationDNF = Vector{Vector{String}}

"""
    Annotations = Dict{String,Vector{String}}

Dictionary used to store (possible multiple) standardized annotations of
something, such as a metabolite or a reaction (as listed by
[`metabolites`](@ref) and [`reactions`](@ref)).

# Example
```
Annotations("PubChem" => ["CID12345", "CID54321"])
```
"""
const Annotations = Dict{String,Vector{String}}

"""
    Notes = Dict{String,Vector{String}}

Free-form notes about something (e.g. a gene listed by [`genes`](@ref)),
categorized by "topic".
"""
const Notes = Dict{String,Vector{String}}
