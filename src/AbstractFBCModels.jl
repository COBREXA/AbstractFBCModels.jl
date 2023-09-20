"""
$(README)
"""
module AbstractFBCModels

using DocStringExtensions
using SparseArrays
using Downloads, SHA, SBML

import RequiredInterfaces as R
import PikaParser as PP

include("types.jl")
include("accessors.jl")
include("io.jl")

# utility functions used to read models, common to many file types
include("utils.jl")

# define required interface for FBCModel types
R.@required AbstractFBCModel begin
    # IO
    # load(::Type{JSONFBCModel}, ::String) # TODO RequiredInterfaces.jl issue #13
    save(::AbstractFBCModel, ::String)
    # filename_extensions(::Type{JSONFBCModel}) # TODO RequiredInterfaces.jl issue #13
    # basic accessors
    reactions(::AbstractFBCModel)
    n_reactions(::AbstractFBCModel)
    metabolites(::AbstractFBCModel)
    n_metabolites(::AbstractFBCModel)
    genes(::AbstractFBCModel)
    n_genes(::AbstractFBCModel)
    # optimization problem
    stoichiometry(::AbstractFBCModel)
    bounds(::AbstractFBCModel)
    balance(::AbstractFBCModel)
    objective(::AbstractFBCModel)
    # reaction annotations
    reaction_gene_associations(::AbstractFBCModel,::String)
    reaction_stoichiometry(::AbstractFBCModel, ::String)
    reaction_annotations(::AbstractFBCModel, ::String)
    reaction_notes(::AbstractFBCModel, ::String)
    reaction_name(::AbstractFBCModel, ::String)
    # metabolite annotations
    metabolite_formula(::AbstractFBCModel, ::String)
    metabolite_charge(::AbstractFBCModel, ::String)
    metabolite_compartment(::AbstractFBCModel, ::String)
    metabolite_annotations(::AbstractFBCModel, ::String)
    metabolite_notes(::AbstractFBCModel, ::String)
    metabolite_name(::AbstractFBCModel, ::String)
    # gene annotations    
    gene_annotations(::AbstractFBCModel, ::String)
    gene_notes(::AbstractFBCModel, ::String)
    gene_name(::AbstractFBCModel, ::String)
    # conversion
    # convert(::Type{JSONFBCModel}, model::AbstractFBCModel) # TODO RequiredInterfaces.jl issue #13
end

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
