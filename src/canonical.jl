
module CanonicalModel

using DocStringExtensions

import ..AbstractFBCModel as A
import Serialization as S

Base.@kwdef mutable struct Reaction
    name::Maybe{String} = nothing
    lower_bound::Float64 = -Inf
    upper_bound::Float64 = Inf
    stoichiometry::Dict{String,Float64} = Dict()
    objective_coefficient::Float64 = 0.0
    gene_association_dnf::Maybe{A.GeneAssociationDNF} = nothing
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

Base.@kwdef mutable struct Metabolite
    name::Maybe{String} = nothing
    compartment::Maybe{String} = nothing
    formula::Maybe{A.MetaboliteFormula} = nothing
    charge::Maybe{Int} = nothing
    balance::Float64 = 0.0
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

Base.@kwdef mutable struct Gene
    name::Maybe{String} = nothing
    annotations::A.Annotations = nothing
    notes::A.Notes = nothing
end

Base.@kwdef struct Model <: AbstractFBCModel
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

A.stoichiometry(m::Model) = sparse(
    get(m.reactions[rid].stoichiometry, mid, 0.0) for mid in A.metabolites(m),
    rid in A.reactions(m)
)

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

function A.load(::Type{Model}, path::String)::Model
    S.deserialize(path)
end

A.save(m::Model, path::String) = S.serialize(path, m)

A.filename_extensions(::Type{Model}) = ["canonical-serialized-fbc"]

end # module CanonicalModel
