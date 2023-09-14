
_missing_impl_error(m, a) = throw(MethodError(m, a))

"""
$(TYPEDSIGNATURES)

Return a vector of reaction identifiers in a model. The vector precisely
corresponds to the columns in [`stoichiometry`](@ref) matrix.

For technical reasons, the "reactions" may sometimes not be true reactions but
various virtual and helper pseudo-reactions that are used in the metabolic
modeling, such as metabolite exchanges, separate forward and reverse reactions,
supplies of enzymatic and genetic material and virtual cell volume, etc. To
simplify the view of the model contents use [`reaction_flux`](@ref).
"""
function reactions(a::MetabolicModel)::Vector{String}
    _missing_impl_error(reactions, (a,))
end

"""
$(TYPEDSIGNATURES)

Return a vector of metabolite identifiers in a model. The vector precisely
corresponds to the rows in [`stoichiometry`](@ref) matrix.

As with [`reactions`](@ref)s, some metabolites in models may be virtual,
representing purely technical equality constraints.
"""
function metabolites(a::MetabolicModel)::Vector{String}
    _missing_impl_error(metabolites, (a,))
end

"""
$(TYPEDSIGNATURES)

Get the number of reactions in a model.
"""
function n_reactions(a::MetabolicModel)::Int
    length(reactions(a))
end

"""
$(TYPEDSIGNATURES)

Get the number of metabolites in a model.
"""
function n_metabolites(a::MetabolicModel)::Int
    length(metabolites(a))
end

"""
$(TYPEDSIGNATURES)

Get the sparse stoichiometry matrix of a model. A feasible solution `x` of a
model `m` is defined as satisfying the equations:

- `stoichiometry(m) * x .== balance(m)`
- `x .>= lbs`
- `y .<= ubs`
- `(lbs, ubs) == bounds(m)
"""
function stoichiometry(a::MetabolicModel)::SparseMat
    _missing_impl_error(stoichiometry, (a,))
end

"""
$(TYPEDSIGNATURES)

Get the lower and upper solution bounds of a model.
"""
function bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}
    _missing_impl_error(bounds, (a,))
end

"""
$(TYPEDSIGNATURES)

Get the sparse balance vector of a model.
"""
function balance(a::MetabolicModel)::SparseVec
    return spzeros(n_metabolites(a))
end

"""
$(TYPEDSIGNATURES)

Get the objective vector of the model. Analysis functions, such as
[`flux_balance_analysis`](@ref), are supposed to maximize `dot(objective, x)`
where `x` is a feasible solution of the model.
"""
function objective(a::MetabolicModel)::SparseVec
    _missing_impl_error(objective, (a,))
end

"""
$(TYPEDSIGNATURES)

Get a matrix of coupling constraint definitions of a model. By default, there
is no coupling in the models.
"""
function coupling(a::MetabolicModel)::SparseMat
    return spzeros(0, n_reactions(a))
end

"""
$(TYPEDSIGNATURES)

Get the number of coupling constraints in a model.
"""
function n_coupling_constraints(a::MetabolicModel)::Int
    size(coupling(a), 1)
end

"""
$(TYPEDSIGNATURES)

Get the lower and upper bounds for each coupling bound in a model, as specified
by `coupling`. By default, the model does not have any coupling bounds.
"""
function coupling_bounds(a::MetabolicModel)::Tuple{Vector{Float64},Vector{Float64}}
    return (spzeros(0), spzeros(0))
end

"""
$(TYPEDSIGNATURES)

Return identifiers of all genes contained in the model. By default, there are
no genes.

In SBML, these are usually called "gene products" but we write `genes` for
simplicity.
"""
function genes(a::MetabolicModel)::Vector{String}
    return []
end

"""
$(TYPEDSIGNATURES)

Return the number of genes in the model (as returned by [`genes`](@ref)). If
you just need the number of the genes, this may be much more efficient than
calling [`genes`](@ref) and measuring the array.
"""
function n_genes(a::MetabolicModel)::Int
    return length(genes(a))
end

"""
$(TYPEDSIGNATURES)

Returns the sets of genes that need to be present so that the reaction can work
(technically, a DNF on gene availability, with positive atoms only).

For simplicity, `nothing` may be returned, meaning that the reaction always
takes place. (in DNF, that would be equivalent to returning `[[]]`.)
"""
function reaction_gene_association(
    a::MetabolicModel,
    reaction_id::String,
)::Maybe{GeneAssociation}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the subsystem of reaction `reaction_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function reaction_subsystem(model::MetabolicModel, reaction_id::String)::Maybe{String}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid` in the model. The dictionary
maps the metabolite IDs to their stoichiometric coefficients.
"""
function reaction_stoichiometry(m::MetabolicModel, rid::String)::Dict{String,Float64}
    mets = metabolites(m)
    Dict(
        mets[k] => v for
        (k, v) in zip(findnz(stoichiometry(m)[:, first(indexin([rid], reactions(m)))])...)
    )
end

"""
$(TYPEDSIGNATURES)

Return the formula of metabolite `metabolite_id` in `model`.
Return `nothing` in case the formula is not known or irrelevant.
"""
function metabolite_formula(
    model::MetabolicModel,
    metabolite_id::String,
)::Maybe{MetaboliteFormula}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the charge associated with metabolite `metabolite_id` in `model`.
Returns `nothing` if charge not present.
"""
function metabolite_charge(model::MetabolicModel, metabolite_id::String)::Maybe{Int}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the compartment of metabolite `metabolite_id` in `model` if it is assigned. If not,
return `nothing`.
"""
function metabolite_compartment(model::MetabolicModel, metabolite_id::String)::Maybe{String}
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the name of reaction with ID `rid`.
"""
reaction_name(model::MetabolicModel, rid::String) = nothing

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `mid`.
"""
metabolite_name(model::MetabolicModel, mid::String) = nothing

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `gid`.
"""
gene_name(model::MetabolicModel, gid::String) = nothing
