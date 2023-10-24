
"""
$(TYPEDSIGNATURES)

Return a vector of reaction identifiers in a model.

For technical reasons, the "reactions" may sometimes not be true reactions but
various virtual and helper pseudo-reactions that are used in the metabolic
modeling, such as metabolite exchanges, separated forward and reverse reaction
directions, supplies of enzymatic and genetic material and virtual cell volume,
etc.
"""
reactions(a::AbstractFBCModel)::Vector{String} = unimplemented(typeof(a), :reactions)

"""
$(TYPEDSIGNATURES)

The number of reactions in the given model. Must be equal to the length of the
vector returned by [`reactions`](@ref), and may be more efficient for just
determining the size.
"""
n_reactions(a::AbstractFBCModel)::Int = unimplemented(typeof(a), :n_reactions)

"""
$(TYPEDSIGNATURES)

Return a vector of metabolite identifiers in a model.

As with [`reactions`](@ref), some metabolites in models may be virtual,
representing purely technical equality constraints.
"""
metabolites(a::AbstractFBCModel)::Vector{String} = unimplemented(typeof(a), :metabolites)

"""
$(TYPEDSIGNATURES)

The number of metabolites in the given model. Must be equal to the length of
the vector returned by [`metabolites`](@ref), and may be more efficient for
just determining the size.
"""
function n_metabolites end
n_metabolites(a::AbstractFBCModel)::Int = unimplemented(typeof(a), :n_metabolites)

"""
$(TYPEDSIGNATURES)

Return identifiers of all genes contained in the model. Empty if none.

Genes are also sometimes called "gene products" but we write genes for
simplicity.
"""
genes(::AbstractFBCModel)::Vector{String} = unimplemented(typeof(a), :genes)

"""
$(TYPEDSIGNATURES)

The number of genes in the model (must be equal to the length of vector given
by [`genes`](@ref)).

This may be more efficient than calling [`genes`](@ref) and measuring the
array.
"""
n_genes(::AbstractFBCModel)::Int = unimplemented(typeof(a), :n_genes)

"""
$(TYPEDSIGNATURES)

The sparse stoichiometric matrix of a given model.

This usually corresponds to all the equality constraints in the model. The
matrix must be of size [`n_metabolites`](@ref) by [`n_reactions`](@ref).
"""
stoichiometry(a::AbstractFBCModel)::SparseMat = unimplemented(typeof(a), :stoichiometry)

"""
$(TYPEDSIGNATURES)

Get the lower and upper bounds of all reactions in a model.
"""
bounds(a::AbstractFBCModel)::Tuple{Vector{Float64},Vector{Float64}} =
    unimplemented(typeof(a), :bounds)

"""
$(TYPEDSIGNATURES)

Get the sparse balance vector of a model, which usually corresponds to the
accumulation term associated with stoichiometric matrix.
"""
balance(a::AbstractFBCModel)::SparseVec = spzeros(length(n_metabolites(a)))

"""
$(TYPEDSIGNATURES)

Get the objective vector of the model.
"""
objective(a::AbstractFBCModel)::SparseVec = unimplemented(typeof(a), :objective)

"""
$(TYPEDSIGNATURES)

Evaluate whether the reaction can work given in a conditions given by the
current availability of gene products, or `nothing` if the information is not
recorded. The availability us queried via `gene_products_available`, which must
be a function of a single `String` argument that returns `Bool`.

Generally it may be simpler to use [`reaction_gene_association_dnf`](@ref), but
in many models the complexity of the conversion to DNF is prohibitive.
"""
reaction_gene_products_available(
    ::AbstractFBCModel,
    reaction_id::String,
    gene_product_available::Function,
)::Maybe{Bool} = nothing

"""
$(TYPEDSIGNATURES)

Returns the sets of genes that need to be present for the reaction to work in a
DNF formula. This helps for constructively using the reaction-gene-association
information.
"""
reaction_gene_association_dnf(
    ::AbstractFBCModel,
    reaction_id::String,
)::Maybe{GeneAssociationDNF} = nothing

"""
$(TYPEDSIGNATURES)

The stoichiometry of the given reaction as a dictionary
maps the metabolite IDs to their stoichiometric coefficients.

Using this function may be more efficient in cases than loading the whole
[`stoichiometry`](@ref).
"""
reaction_stoichiometry(a::AbstractFBCModel, reaction_id::String)::Dict{String,Float64} =
    unimplemented(typeof(a), :reaction_stoichiometry)

"""
$(TYPEDSIGNATURES)

A dictionary of standardized names that may help identifying the reaction. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. `"Reactome" => ["reactomeID123"]`.
"""
reaction_annotations(::AbstractFBCModel, reaction_id::String)::Annotations = Dict()

"""
$(TYPEDSIGNATURES)

Free-text notes organized in a dictionary by topics about the given reaction in
the model.
"""
reaction_notes(::AbstractFBCModel, reaction_id::String)::Notes = Dict()

"""
$(TYPEDSIGNATURES)

Name of the given reaction.
"""
reaction_name(::AbstractFBCModel, reaction_id::String)::Maybe{String} = nothing


"""
$(TYPEDSIGNATURES)

The formula of the given metabolite in the model, or `nothing` in case the
formula is not recorded.
"""
metabolite_formula(::AbstractFBCModel, metabolite_id::String)::Maybe{MetaboliteFormula} =
    nothing

"""
$(TYPEDSIGNATURES)

The charge of the given metabolite in the model, or `nothing` in case the
charge is not recorded.
"""
metabolite_charge(::AbstractFBCModel, metabolite_id::String)::Maybe{Int} = nothing

"""
$(TYPEDSIGNATURES)

The compartment of the given metabolite in the model. `nothing` if no
compartment is assigned.
"""
metabolite_compartment(::AbstractFBCModel, metabolite_id::String)::Maybe{String} = nothing

"""
$(TYPEDSIGNATURES)

A dictionary of standardized names that may help to identify the metabolite. The
dictionary should assigns vectors of possible identifiers to identifier system
names, such as `"ChEMBL" => ["123"]` or
`"PubChem" => ["CID123", "CID654645645"]`.
"""
metabolite_annotations(::AbstractFBCModel, metabolite_id::String)::Annotations = Dict()

"""
$(TYPEDSIGNATURES)

Free-text notes organized in a dictionary by topics about the given metabolite
in the model.
"""
metabolite_notes(::AbstractFBCModel, metabolite_id::String)::Notes = Dict()

"""
$(TYPEDSIGNATURES)

The name of the given metabolite, if assigned.
"""
metabolite_name(::AbstractFBCModel, metabolite_id::String)::Maybe{String} = nothing

"""
$(TYPEDSIGNATURES)

A dictionary of standardized names that identify the corresponding gene or
product. The dictionary assigns vectors of possible identifiers to identifier
system names, such as `"PDB" => ["PROT01"]`.
"""
gene_annotations(::AbstractFBCModel, gene_id::String)::Annotations = Dict()

"""
$(TYPEDSIGNATURES)

Free-text notes organized in a dictionary by topics about the given gene in the
model.
"""
gene_notes(::AbstractFBCModel, gene_id::String)::Notes = Dict()

"""
$(TYPEDSIGNATURES)

The name of the given gene in the model, if recorded.
"""
gene_name(::AbstractFBCModel, gene_id::String)::Maybe{String} = nothing
