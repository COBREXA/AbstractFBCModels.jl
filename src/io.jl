
"""
    load(a::Type{AbstractFBCModel}, path::String)::Nothing

Load a flux balance model from path.
"""
function load end

"""
    save(a::AbstractFBCModel, path::String)::Nothing

Save a flux balance model to path.
"""
function save end

"""
    filename_extensions(a::Type{AbstractFBCModel})::Vector{String}

This function declares a vector of filename extensions that are common for files
that contain the given metabolic model. This is used by [load](@ref) to guess
the type of the model to load.
"""
function filename_extensions end

# """
# $(TYPEDSIGNATURES)
# """
# function guess_model_type_from_filename(path::String)
#     # TODO
#     error("filename extension not recognized")
# end

# """
# $(TYPEDSIGNATURES)

# Load a model from path. The type of the model is automatically guessed based
# on the filename extension.
# """
# load(path::String) = load(guess_model_type_from_filename(path), path)

Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractFBCModel) =
    print(
        io,
        "$(typeof(m))(#= $(n_reactions(m)) reactions, $(n_metabolites(m)) metabolites =#)",
    )
