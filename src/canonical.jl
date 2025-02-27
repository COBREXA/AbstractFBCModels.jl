
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

A canonical Julia representation of a row in a coupling matrix of the
`AbstractFBCModels` interface.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Coupling
    name::A.Maybe{String} = nothing
    reaction_weights::Dict{String,Float64} = Dict()
    lower_bound::Float64 = -Inf
    upper_bound::Float64 = Inf
    annotations::A.Annotations = A.Annotations()
    notes::A.Notes = A.Notes()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Coupling) = A.pretty_print_kwdef(io, x)

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
Base.@kwdef mutable struct Model <: A.AbstractFBCModel
    reactions::Dict{String,Reaction} = Dict()
    metabolites::Dict{String,Metabolite} = Dict()
    genes::Dict{String,Gene} = Dict()
    couplings::Dict{String,Coupling} = Dict()
end

Base.show(io::Base.IO, ::MIME"text/plain", x::Model) = A.pretty_print_kwdef(io, x)

A.reactions(m::Model) = sort(collect(keys(m.reactions)))
A.metabolites(m::Model) = sort(collect(keys(m.metabolites)))
A.genes(m::Model) = sort(collect(keys(m.genes)))
A.couplings(m::Model) = sort(collect(keys(m.couplings)))
A.n_reactions(m::Model) = length(m.reactions)
A.n_metabolites(m::Model) = length(m.metabolites)
A.n_genes(m::Model) = length(m.genes)
A.n_couplings(m::Model) = length(m.couplings)
A.reaction_name(m::Model, id::String) = m.reactions[id].name
A.metabolite_name(m::Model, id::String) = m.metabolites[id].name
A.gene_name(m::Model, id::String) = m.genes[id].name
A.coupling_name(m::Model, id::String) = m.couplings[id].name
A.reaction_annotations(m::Model, id::String) = m.reactions[id].annotations
A.metabolite_annotations(m::Model, id::String) = m.metabolites[id].annotations
A.gene_annotations(m::Model, id::String) = m.genes[id].annotations
A.coupling_annotations(m::Model, id::String) = m.couplings[id].annotations
A.reaction_notes(m::Model, id::String) = m.reactions[id].notes
A.metabolite_notes(m::Model, id::String) = m.metabolites[id].notes
A.gene_notes(m::Model, id::String) = m.genes[id].notes
A.coupling_notes(m::Model, id::String) = m.couplings[id].notes

function A.stoichiometry(m::Model)
    midxs = Dict(mid => idx for (idx, mid) in enumerate(A.metabolites(m)))
    I = Int[]
    J = Int[]
    V = Float64[]
    for (ridx, rid) in enumerate(A.reactions(m))
        for (smid, v) in m.reactions[rid].stoichiometry
            push!(I, midxs[smid])
            push!(J, ridx)
            push!(V, v)
        end
    end
    sparse(I, J, V, A.n_metabolites(m), A.n_reactions(m))
end

function A.coupling(m::Model)
    ridxs = Dict(rid => idx for (idx, rid) in enumerate(A.reactions(m)))
    I = Int[]
    J = Int[]
    V = Float64[]
    for (cidx, cid) in enumerate(A.couplings(m))
        for (rid, v) in m.couplings[cid].reaction_weights
            push!(I, cidx)
            push!(J, ridxs[rid])
            push!(V, v)
        end
    end
    sparse(I, J, V, A.n_couplings(m), A.n_reactions(m))
end

A.bounds(m::Model) = (
    [m.reactions[rid].lower_bound for rid in A.reactions(m)],
    [m.reactions[rid].upper_bound for rid in A.reactions(m)],
)

A.coupling_bounds(m::Model) = (
    [m.couplings[cid].lower_bound for cid in A.couplings(m)],
    [m.couplings[cid].upper_bound for cid in A.couplings(m)],
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

A.coupling_weights(m::Model, id::String) = m.couplings[id].reaction_weights

A.load(::Type{Model}, path::String) = S.deserialize(path)
A.save(m::Model, path::String) = S.serialize(path, m)
A.filename_extensions(::Type{Model}) = ["canonical-serialized-fbc"]

function Base.convert(::Type{Model}, x::A.AbstractFBCModel)
    (lbs, ubs) = A.bounds(x)
    (clbs, cubs) = A.coupling_bounds(x)
    Model(
        reactions = Dict(
            r => Reaction(
                name = A.reaction_name(x, r),
                lower_bound = lb,
                upper_bound = ub,
                stoichiometry = A.reaction_stoichiometry(x, r),
                objective_coefficient = o,
                gene_association_dnf = A.reaction_gene_association_dnf(x, r),
                annotations = A.reaction_annotations(x, r),
                notes = A.reaction_notes(x, r),
            ) for (r, o, lb, ub) in zip(A.reactions(x), A.objective(x), lbs, ubs)
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
            ) for (m, b) in zip(A.metabolites(x), A.balance(x))
        ),
        genes = Dict(
            g => Gene(
                name = A.gene_name(x, g),
                annotations = A.gene_annotations(x, g),
                notes = A.gene_notes(x, g),
            ) for g in A.genes(x)
        ),
        couplings = Dict(
            c => Coupling(
                name = A.coupling_name(x, c),
                lower_bound = lb,
                upper_bound = ub,
                reaction_weights = A.coupling_weights(x, c),
                annotations = A.coupling_annotations(x, c),
                notes = A.coupling_notes(x, c),
            ) for (c, lb, ub) in zip(A.couplings(x), clbs, cubs)
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Rebuild the given [`Model`](@ref) with all identifiers changed accordingly to
the given individual identifier-mapping functions in arguments. These take a
single `String` identifier as an argument, and return a new `String`
identifier. By default, the identifier maps are `identity`, i.e., no
identifiers are changed.

The identifier maps *must be injective and pure within `identifier_map`*, i.e.,
they must not create identifier name conflicts, and the results returned for a
given argument must be the same during the whole call of `identifier_map`.
Errors stemming from use of non-injective and impure identifier maps are not
handled.

Internal data structures are copied only by reference wherever possible.
"""
function identifier_map(
    m::Model;
    reaction_map = identity,
    metabolite_map = identity,
    gene_map = identity,
    coupling_map = identity,
    compartment_map = identity,
)
    return Model(;
        reactions = Dict(
            reaction_map(k) => Reaction(;
                name = v.name,
                lower_bound = v.lower_bound,
                upper_bound = v.upper_bound,
                stoichiometry = Dict(
                    metabolite_map(sk) => sv for (sk, sv) in v.stoichiometry
                ),
                objective_coefficient = v.objective_coefficient,
                gene_association_dnf = [
                    String[gene_map(g) for g in gs] for gs in a.gene_association_dnf
                ],
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.reactions
        ),
        metabolites = Dict(
            metabolite_map(k) => Metabolite(;
                name = v.name,
                compartment = isnothing(v.compartment) ? nothing :
                              compartment_map(v.compartment),
                formula = v.formula,
                charge = v.charge,
                balance = v.balance,
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.metabolites
        ),
        genes = Dict(gene_map(k) => v for (k, v) in m.genes),
        couplings = Dict(
            coupling_map(k) => Coupling(;
                name = v.name,
                reaction_weights = Dict(
                    reaction_map(rk) => rv for (rk, rv) in v.reaction_weights
                ),
                lower_bound = v.lower_bound,
                upper_bound = v.upper_bound,
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.couplings
        ),
    )
end

end # module CanonicalModel
