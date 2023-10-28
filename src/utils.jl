
import Downloads: download
import SHA: sha256
import InteractiveUtils: methodswith

unimplemented(t::Type, x::Symbol) =
    error("AbstractFBCModels interface method $x is not implemented for type $t")

"""
$(TYPEDSIGNATURES)

Provide a `methodswith`-style listing of accessors that the model implementors
should implement.

For typesystem reasons, the list **will not contain** methods for
[`save`](@ref) and [`filename_extensions`](@ref) that dispatch on type objects.
You should implement these as well.
"""
function accessors()
    ms = Method[]
    for nm in names(AbstractFBCModels; all = true)
        f = getfield(AbstractFBCModels, nm)
        if isa(f, Base.Callable)
            try
                methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
            catch
                # Earlier versions of Julia tend to throw MethodErrors here
                # whenever the method actually doesn't exist (e.g. 1.6.x
                # reports that it's actually a missing `methodswith` method
                # rather than the one of `f`. If that happens, we can simply do
                # nothing.
            end
        end
    end
    ms
end

"""
$(TYPEDSIGNATURES)

Check if the file at `path` has the expected hash.
"""
function check_cached_file_hash(path, expected_checksum)
    actual_checksum = bytes2hex(sha256(open(path)))
    if actual_checksum != expected_checksum
        @warn "The downloaded data file `$path' seems to be different from the expected one. Tests will likely fail." actual_checksum expected_checksum
        @info "You can likely remove `$path' to force re-downloading."
    end
end

"""
$(TYPEDSIGNATURES)

Download the file at `url` and save it at `path`, also check if this file is
the expected file by calling [`check_cached_file_hash`](@ref). If the file has
already been downloaded, and stored at `path`, then it is not downloaded again.
"""
function download_data_file(url, path, hash)
    if isfile(path)
        check_cached_file_hash(path, hash)
        @info "using cached `$path'"
        return path
    end

    path = download(url, path)
    check_cached_file_hash(path, hash)
    return path
end

"""
$(TYPEDSIGNATURES)

Utility function to compute the value of
[`reaction_gene_products_available`](@ref) for models that already implement
[`reaction_gene_association_dnf`](@ref).
"""
function reaction_gene_products_available_from_dnf(
    m::AbstractFBCModel,
    reaction_id::String,
    available::Function,
)::Maybe{Bool}
    gss = reaction_gene_association_dnf(m, reaction_id)
    isnothing(gss) ? nothing : any(gs -> all(available, gs), gss)
end
