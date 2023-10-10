
module CanonicalModel

using DocStringExtensions

import ..AbstractFBCModels as A
import Serialization as S
import SparseArrays: sparse

"""
$(TYPEDEF)

A canonical Julia representation of a reaction in the `AbstractFBCModels` interface.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Reaction
    name::A.Maybe{String} = nothing
    lower_bound::Float64 = -Inf
    upper_bound::Float64 = Inf
    stoichiometry::Dict{String,Float64} = Dict()
    objective_coefficient::Float64 = 0.0
    gene_association_dnf::A.Maybe{A.GeneAssociationDNF} = nothing
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

"""
$(TYPEDEF)

A canonical Julia representation of a metabolite in the `AbstractFBCModels` interface.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Metabolite
    name::A.Maybe{String} = nothing
    compartment::A.Maybe{String} = nothing
    formula::A.Maybe{A.MetaboliteFormula} = nothing
    charge::A.Maybe{Int} = nothing
    balance::Float64 = 0.0
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

"""
$(TYPEDEF)

A canonical Julia representation of a gene in the `AbstractFBCModels` interface.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Gene
    name::A.Maybe{String} = nothing
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

"""
$(TYPEDEF)

A canonical Julia representation of a metabolic model that sotres exactly the
data represented by `AbstractFBCModels` accessors.

The implementation is useful for manipulating model data manually without
writing new model types, or even for constructing models from base data in many
simple cases.

Additionally, you can use the implementation of accessors for this model type
in the source code of `AbstractFBCModels` as a starting point for creating an
`AbstractFBCModel` interface for your own models.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct Model <: A.AbstractFBCModel
    reactions::Dict{String,Reaction}
    metabolites::Dict{String,Metabolite}
    genes::Dict{String,Gene}
end

A.reactions(m::Model) = sort(collect(keys(m.reactions)))
A.metabolites(m::Model) = sort(collect(keys(m.metabolites)))
A.genes(m::Model) = sort(collect(keys(m.metabolites)))
A.n_reactions(m::Model) = length(m.reactions)
A.n_metabolites(m::Model) = length(m.metabolites)
A.n_genes(m::Model) = length(m.genes)
A.reaction_name(m::Model, id::String) = m.reactions[id].name
A.metabolite_name(m::Model, id::String) = m.metabolites[id].name
A.gene_name(m::Model, id::String) = m.genes[id].name
A.reaction_annotations(m::Model, id::String) = m.reactions[id].annotations
A.metabolite_annotations(m::Model, id::String) = m.metabolites[id].annotations
A.gene_annotations(m::Model, id::String) = m.genes[id].annotations
A.reaction_notes(m::Model, id::String) = m.reactions[id].notes
A.metabolite_notes(m::Model, id::String) = m.metabolites[id].notes
A.gene_notes(m::Model, id::String) = m.genes[id].notes

A.stoichiometry(m::Model) =
    let rids = A.reactions(m)
        sparse(
            get(m.reactions[rid].stoichiometry, mid, 0.0) for mid in A.metabolites(m),
            rid in rids
        )
    end

A.bounds(m::Model) = (
    [m.reactions[rid].lower_bound for rid in A.reactions(m)],
    [m.reactions[rid].upper_bound for rid in A.reactions(m)],
)

A.balance(m::Model) = sparse(m.metabolite[mid].balance for mid in A.metabolites(m))
A.objective(m::Model) =
    sparse(m.reaction[rid].objective_coefficient for rid in A.reactions(m))

A.reaction_gene_association_dnf(m::Model, id::String) = m.reactions[id].gene_association_dnf
A.reaction_gene_products_available(m::Model, id::String, fn::Function) =
    A.reaction_gene_products_available_from_dnf(m, id, fn)
A.reaction_stoichiometry(m::Model, id::String) = m.reactions[id].stoichiometry

A.metabolite_formula(m::Model, id::String) = m.metabolites[id].formula
A.metabolite_charge(m::Model, id::String) = m.metabolites[id].charge
A.metabolite_compartment(m::Model, id::String) = m.metabolites[id].compartment

A.load(::Type{Model}, path::String)::Model = S.deserialize(path)
A.save(m::Model, path::String) = S.serialize(path, m)
A.filename_extensions(::Type{Model}) = ["canonical-serialized-fbc"]

end # module CanonicalModel
