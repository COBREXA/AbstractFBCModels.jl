
"""
    reactions(a::AbstractFBCModel)::Vector{String}

Return a vector of reaction identifiers in a model.

For technical reasons, the "reactions" may sometimes not be true reactions but
various virtual and helper pseudo-reactions that are used in the metabolic
modeling, such as metabolite exchanges, separate forward and reverse reactions,
supplies of enzymatic and genetic material and virtual cell volume, etc.
"""
function reactions end

"""
    n_reactions(a::AbstractFBCModel)::Int

Get the number of reactions in a model.
"""
function n_reactions end

"""
    metabolites(a::AbstractFBCModel)::Vector{String}

Return a vector of metabolite identifiers in a model. 

As with [reactions](@ref)s, some metabolites in models may be virtual,
representing purely technical equality constraints.
"""
function metabolites end

"""
    n_metabolites(a::AbstractFBCModel)::Int

Get the number of metabolites in a model.
"""
function n_metabolites end

"""
    genes(a::AbstractFBCModel)::Vector{String}

Return identifiers of all genes contained in the model. Empty if none.

Genes are also sometimes called "gene products" but we write genes for
simplicity.
"""
function genes end

"""
    n_genes(a::AbstractFBCModel)::Int

Return the number of genes in the model (as returned by [genes](@ref)). If
you just need the number of the genes, this may be much more efficient than
calling [genes](@ref) and measuring the array.
"""
function n_genes end

"""
    stoichiometry(a::AbstractFBCModel)::SparseMat

Get the sparse stoichiometric matrix of a model. 

This usually corresponds to all the equality constraints, and has dimensions of
[n_metabolites](@ref) by [n_reactions](@ref).
"""
function stoichiometry end

"""
    bounds(a::AbstractFBCModel)::Tuple{Vector{Float64},Vector{Float64}}

Get the lower and upper bounds of the reactions in a model.
"""
function bounds end

"""
    balance(a::AbstractFBCModel)::SparseVec

Get the sparse balance vector of a model, which usually corresponds to the
accumulation term associated with stoichiometric matrix.
"""
function balance end

"""
    objective(a::AbstractFBCModel)::SparseVec

Get the objective vector of the model.
"""
function objective end

"""
    reaction_gene_associations(a::AbstractFBCModel,reaction_id::String)::Maybe{GeneAssociation}

Returns the sets of genes that need to be present so that the reaction can work
(technically, a DNF on gene availability, with positive atoms only).
"""
function reaction_gene_associations end

"""
    reaction_stoichiometry(m::AbstractFBCModel, rid::String)::Dict{String,Float64}

Return the stoichiometry of reaction with ID rid in the model. The dictionary
maps the metabolite IDs to their stoichiometric coefficients.
"""
function reaction_stoichiometry end

"""
    reaction_annotations(a::AbstractFBCModel, reaction_id::String)::Annotations

Return standardized names that may help identifying the reaction. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. "Reactome" => ["reactomeID123"].
"""
function reaction_annotations end

"""
    reaction_notes(model::AbstractFBCModel, reaction_id::String)::Notes

Return the notes associated with reaction reaction_id in model.
"""
function reaction_notes end

"""
    reaction_name(model::AbstractFBCModel, rid::String)::Maybe{String}

Return the name of reaction with ID rid.
"""
function reaction_name end


"""
    metabolite_formula(model::AbstractFBCModel, metabolite_id::String)::Maybe{MetaboliteFormula}

Return the formula of metabolite metabolite_id in model.
Return nothing in case the formula is not known or irrelevant.
"""
function metabolite_formula end

"""
    metabolite_charge(model::AbstractFBCModel, metabolite_id::String)::Maybe{Int}

Return the charge associated with metabolite metabolite_id in model.
Returns nothing if charge not present.
"""
function metabolite_charge end

"""
    metabolite_compartment(model::AbstractFBCModel,metabolite_id::String)::Maybe{String}

Return the compartment of metabolite metabolite_id in model if it is assigned. If not,
return nothing.
"""
function metabolite_compartment end

"""
    metabolite_annotations(a::AbstractFBCModel, metabolite_id::String,)::Annotations

Return standardized names that may help to reliably identify the metabolite. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. "ChEMBL" => ["123"] or "PubChem" => ["CID123", "CID654645645"].
"""
function metabolite_annotations end

"""
    metabolite_notes(model::AbstractFBCModel, metabolite_id::String)::Notes

Return the notes associated with metabolite reaction_id in model.
"""
function metabolite_notes end

"""
    metabolite_name(model::AbstractFBCModel, mid::String)::Maybe{String}

Return the name of metabolite with ID mid.
"""
function metabolite_name end

"""
    gene_annotations(a::AbstractFBCModel, gene_id::String)::Annotations

Return standardized names that identify the corresponding gene or product. The
dictionary assigns vectors of possible identifiers to identifier system names,
e.g. "PDB" => ["PROT01"].
"""
function gene_annotations end

"""
    gene_notes(model::AbstractFBCModel, gene_id::String)::Notes

Return the notes associated with the gene gene_id in model.
"""
function gene_notes end

"""
    gene_name(model::AbstractFBCModel, gid::String)::Maybe{String}

Return the name of gene with ID gid.
"""
function gene_name end

"""
    convert(T::Type{AbstractFBCModel}, model::AbstractFBCModel)::M where M <: AbstractFBCModel

Convert model to type T.
"""
function convert end

