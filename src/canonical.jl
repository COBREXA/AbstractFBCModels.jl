
module CanonicalModel

using DocStringExtensions

import ..AbstractFBCModels as A
import Serialization as S
import SparseArrays: sparse, findnz

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
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Reaction) = A.pretty_print_kwdef(io, x)

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
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Metabolite) = A.pretty_print_kwdef(io, x)

"""
$(TYPEDEF)

A canonical Julia representation of a gene in the `AbstractFBCModels` interface.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Gene
    name::A.Maybe{String} = nothing
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Gene) = A.pretty_print_kwdef(io, x)

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
    reactions::Dict{String,Reaction} = Dict()
    metabolites::Dict{String,Metabolite} = Dict()
    genes::Dict{String,Gene} = Dict()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Model) = A.pretty_print_kwdef(io, x)

A.reactions(m::Model) = sort(collect(keys(m.reactions)))
A.metabolites(m::Model) = sort(collect(keys(m.metabolites)))
A.genes(m::Model) = sort(collect(keys(m.genes)))
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

function A.stoichiometry(m::Model)
    midxs = Dict(mid => idx for (idx,(mid,_)) = enumerate(m.metabolites))
    I = Int[]
    J = Int[]
    V = Float64[]
    for (ridx, (_, r)) = enumerate(m.reactions)
        for (smid, v) = r.stoichiometry
            push!(I, midxs[smid])
            push!(J, ridx)
            push!(V, v)
        end
    end
    sparse(I, J, V, length(m.metabolites), length(m.reactions))
end

A.bounds(m::Model) = (
    [m.reactions[rid].lower_bound for rid in A.reactions(m)],
    [m.reactions[rid].upper_bound for rid in A.reactions(m)],
)

A.balance(m::Model) =
    sparse(Float64[m.metabolites[mid].balance for mid in A.metabolites(m)])
A.objective(m::Model) =
    sparse(Float64[m.reactions[rid].objective_coefficient for rid in A.reactions(m)])

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

function Base.convert(::Type{Model}, x::A.AbstractFBCModel)
    (lbs, ubs) = A.bounds(x)
    os = A.objective(x)
    bs = A.balance(x)
    mets = A.metabolites(x)
    Model(
        reactions = Dict(
            r => Reaction(
                name = A.reaction_name(x, r),
                lower_bound = lb,
                upper_bound = ub,
                stoichiometry = A.reaction_stoichiometry(x, r),
                gene_association_dnf = A.reaction_gene_association_dnf(x, r),
                annotations = A.reaction_annotations(x, r),
                notes = A.reaction_notes(x, r),
            ) for (r, lb, ub) in zip(A.reactions(x), lbs, ubs)
        ),
        metabolites = Dict(
            m => Metabolite(
                name = A.metabolite_name(x, m),
                balance = b,
                formula = A.metabolite_formula(x, m),
                charge = A.metabolite_charge(x, m),
                compartment = A.metabolite_compartment(x, m),
                annotations = A.metabolite_annotations(x, m),
                notes = A.metabolite_notes(x, m),
            ) for (m, b) in zip(mets, bs)
        ),
        genes = Dict(
            g => Gene(
                name = A.gene_name(x, g),
                annotations = A.gene_annotations(x, g),
                notes = A.gene_notes(x, g),
            ) for g in A.genes(x)
        ),
    )
end

end # module CanonicalModel
