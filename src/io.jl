
"""
$(TYPEDSIGNATURES)

Load a model from path.
"""
function load(::Type{A}, path::String) where A <: AbstractFBCModel
    unimplemented(A, :load)
end

"""
$(TYPEDSIGNATURES)

Save a model to the given path.
"""
save(a::AbstractFBCModel, path::String)::Nothing = unimplemented(typeof(a), :save)

"""
$(TYPEDSIGNATURES)

A vector of filename extensions that are common for files that contain the
given metabolic model type. This is used by [`load`](@ref) to guess the type of
the model that should be loaded.
"""
function filename_extensions(::Type{A}) where A <: AbstractFBCModel
    unimplemented(A, :filename_extensions)
end

"""
$(TYPEDSIGNATURES)

Find which of the known subtypes of [`AbstractFBCModel`](@ref) would typically
be able to open the file at `path`, using information in `filename_extensions`.
"""
eligible_model_types_for_filename(path::String) =
    [t for t = subtypes(A) if any(x -> endswith(".$x", path), filename_extensions(t))]

"""
$(TYPEDSIGNATURES)

Guess which of the known subtypes of [`AbstractFBCModel`](@ref) would typically
open the path (internally using [`eligible_model_types_for_filename`](@ref) for
the purpose). Throws an error if the match is ambiguous or missing.
"""
function guess_model_type_from_filename(path::String)
    ts = eligible_model_types_for_filename(path)
    if length(ts) == 1
        return first(ts)
    elseif isempty(ts)
        error("filename extension not recognized")
    else
        error("multiple model types match the extension: $ts")
    end
end

"""
$(TYPEDSIGNATURES)

Load a model from path. The type of the model is automatically guessed based on
the filename extension. The guessing inspects all subtypes of
[`AbstractFBCModel`](@ref), thus always requires compilation -- if possible,
specify the proper type using the 2-parameter version of [`load`](@ref) to save
time.
"""
load(path::String) = load(guess_model_type_from_filename(path), path)

Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractFBCModel) = print(
    io,
    "$(typeof(m))(#= $(n_reactions(m)) reactions, $(n_metabolites(m)) metabolites =#)",
)
