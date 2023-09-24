"""
$(README)
"""
module AbstractFBCModels

using DocStringExtensions

include("types.jl")
include("accessors.jl")
include("io.jl")
include("utils.jl")

export load,
    save,
    filename_extensions,
    reactions,
    n_reactions,
    metabolites,
    n_metabolites,
    genes,
    n_genes,
    stoichiometry,
    bounds,
    balance,
    objective,
    reaction_gene_associations,
    reaction_stoichiometry,
    reaction_annotations,
    reaction_notes,
    reaction_name,
    metabolite_formula,
    metabolite_charge,
    metabolite_compartment,
    metabolite_annotations,
    metabolite_notes,
    metabolite_name,
    gene_annotations,
    gene_notes,
    gene_name

end # module AbstractFBCModel
