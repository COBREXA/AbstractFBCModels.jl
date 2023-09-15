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
end

end # module AbstractFBCModel
