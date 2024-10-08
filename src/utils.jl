
import Downloads: download
import SHA: sha256
import InteractiveUtils: methodswith


const REQUIRED_FUNCTIONS = Function[]
macro required(sig)
    call_ex = sig.head == :(::) ? sig.args[1] : sig
    call_ex.head == :call || error("malformed signature definition")
    name = call_ex.args[1]
    name isa Symbol || error("malformed signature definition")
    
    model_arg_ex = call_ex.args[2]
    model_arg_ex.head == :(::) || error("malformed signature definition")
    model_arg = model_arg_ex.args[1]
    model_arg isa Symbol || error("malformed signature definition")

    return esc(quote
        Base.@__doc__ $call_ex = $unimplemented(typeof($model_arg), $(Meta.quot(name))) 
        push!(REQUIRED_FUNCTIONS, $name)
        $name
    end)
end

unimplemented(t::Type, x::Symbol) =
    error("AbstractFBCModels interface method $x is not implemented for type $t")

"""
$(TYPEDSIGNATURES)

Provide a `methodswith`-style listing of accessors that the model implementors
may implement.

For typesystem reasons, the list **will not contain** methods for
[`save`](@ref) and [`filename_extensions`](@ref) that dispatch on type objects.
You should implement these as well.

See also [`required_functions`](@ref) for the minimal list that must be implemented.
"""
function accessors()
    ms = Method[]
    for nm in names(AbstractFBCModels; all = true)
        f = getfield(AbstractFBCModels, nm)
        if isa(f, Base.Callable)
            methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
        end
    end
    return ms
end

"""
$(TYPEDSIGNATURES)

Provide a `methodswith`-style listing of functions that the model implementors
must implement to have a functional `AbstractFBCModel`.
constrast this to the longer list of items that are returned by [`accessors`](@ref).
The extra elements have sensible defaults based on these required functions.
Where-as not defining these methods will result in errors.
Though depending on your models capability relying on those defaults may mean some
functionality is hidden.
(e.g. default [`coupling`](@ref) if you don't implement that is to assume none)

"""
function required_functions()
    ms = Method[]
    for f = REQUIRED_FUNCTIONS
        methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
    end
    return ms
end




"""
$(TYPEDSIGNATURES)

Check if the file at `path` has the expected hash.

At this point, the hash is always SHA256.
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
already been downloaded and stored at `path`, then it is not downloaded again.
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
