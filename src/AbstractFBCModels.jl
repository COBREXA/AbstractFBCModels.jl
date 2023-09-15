"""
$(README)
"""
module AbstractFBCModels

using DocStringExtensions

import RequiredInterfaces as R

include("types.jl")
include("accessors.jl")
include("io.jl")

# define required interface for FBCModel types
R.@required AbstractFBCModel begin
    # IO
    load(::AbstractFBCModel, ::String)
    save(::AbstractFBCModel, ::String)
    filename_extensions(::AbstractFBCModel)
    # model
    reactions(::AbstractFBCModel)
    n_reactions(::AbstractFBCModel)
    metabolites(::AbstractFBCModel)
    n_metabolites(::AbstractFBCModel)
    genes(::AbstractFBCModel)
    n_genes(::AbstractFBCModel)
    stoichiometry(::AbstractFBCModel)
    bounds(::AbstractFBCModel)
    balance(::AbstractFBCModel)
    objective(::AbstractFBCModel)
    reaction_gene_association(::AbstractFBCModel,::String)
    reaction_stoichiometry(::AbstractFBCModel, ::String)
    metabolite_formula(::AbstractFBCModel, ::String)
    metabolite_charge(::AbstractFBCModel, ::String)
end

end # module AbstractFBCModel
