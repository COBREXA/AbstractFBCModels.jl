
import Downloads: download
import SHA: sha256
import InteractiveUtils: methodswith


const REQUIRED_ACCESSORS = Function[]
macro required(sig)
    call_ex = sig.head == :(::) ? sig.args[1] : sig
    call_ex.head == :call || error("malformed signature definition")
    name = call_ex.args[1]
    name isa Symbol || error("malformed signature definition")

    model_arg_ex = call_ex.args[2]
    model_arg_ex.head == :(::) || error("malformed signature definition")
    model_arg = model_arg_ex.args[1]
    model_arg isa Symbol || error("malformed signature definition")

    return esc(
        quote
            Base.@__doc__ $call_ex = $unimplemented(typeof($model_arg), $(Meta.quot(name)))
            push!(REQUIRED_ACCESSORS, $name)
            $name
        end,
    )
end

unimplemented(t::Type, x::Symbol) =
    error("AbstractFBCModels interface method $x is not implemented for type $t")


"""
$(TYPEDSIGNATURES)

Provide a `methodswith`-style listing of accessors that the model implementors
may implement.
"""
function accessors()
    ms = Method[]
    for nm in names(AbstractFBCModels; all = true)
        f = getfield(AbstractFBCModels, nm)
        if isa(f, Base.Callable)
            methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
        end
    end

    append!(ms, _type_accessors())
    return ms
end

function _type_accessors()
    # special case: some "accessors" take the Type argument instead of actual instance
    ms = Method[]
    for f in (AbstractFBCModels.load, AbstractFBCModels.filename_extensions)
        for m in methods(f)
            m.sig isa UnionAll || continue
            # Deep magic: basically this matches on `f(::Type{A},...) where A<:AbstractFBCModel`
            type_param = Base.unwrap_unionall(m.sig).parameters[2].parameters[1].ub
            if type_param == AbstractFBCModels.AbstractFBCModel
                push!(ms, m)
            end
        end
    end

    return ms
end


"""
$(TYPEDSIGNATURES)

Provide a `methodswith`-style listing of functions that the model implementors
must implement to have a functional `AbstractFBCModel`.

The output listing is a subset of the longer list returned by
[`accessors`](@ref). Some of the accessors may have sensible defaults -- e.g.,
it is safe to assume an empty default [`coupling`](@ref) if the functions are
not defined. The accessors which do not possess a natural default and thus
*must be defined* (and trying to use a model without them will almost certainly
cause runtime errors) are exactly the ones listed by this function.

"""
function required_accessors()
    ms = Method[]
    for f in REQUIRED_ACCESSORS
        methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
    end
    append!(ms, _type_accessors())
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
