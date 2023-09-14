
"""
$(TYPEDSIGNATURES)
"""
function load(a::Type{AbstractFBCModel}, path::String)
    _missing_impl_error(load, (a,))
end

"""
$(TYPEDSIGNATURES)
"""
function save(a::AbstractFBCModel, path::String)
    _missing_impl_error(save, (a,))
end

"""
$(TYPEDSIGNATURES)

This function declares a vector of filename extensions that are common for
files that contain the given metabolic model. This is used by [`load`](@ref)
and [`save`](@ref) to guess the type of the model to load.
"""
function filename_extensions(a::Type{AbstractFBCModel})::Vector{String}
    _missing_impl_error(filename_extensions, (a,))
end

"""
$(TYPEDSIGNATURES)
"""
function guess_model_type_from_filename(path::String)
    #TODO
    error("filename extension not recognized")
end

"""
$(TYPEDSIGNATURES)

Load a model from `path`. The type of the model is automatically guessed based
on the filename extension.
"""
load(path::String) = load(guess_model_type_from_filename(path), path)
