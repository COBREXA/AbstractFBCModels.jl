
import Downloads: download
import SHA: sha256
import Test: @testset, @test
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
            methodswith(AbstractFBCModels.AbstractFBCModel, f, ms)
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
        @error "The downloaded data file `$path' seems to be different from the expected one. Tests will likely fail." actual_checksum expected_checksum
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

    Downloads.download(url, path)
    check_cached_file_hash(path, hash)
    return path
end

"""
$(TYPEDSIGNATURES)

Test if the given model type loads properly from a file.

The function uses the testing infrastructure from `Test` to report problems --
it is supposed to be a part of larger test-sets, preferably in all model
implementation packages.
"""
function run_fbcmodel_tests(::Type{X}, path::String) where {X<:AbstractFBCModel}
    @testset "Model type $X in file $path" begin
        # TODO optionally download the file as given by kwargs

        model = load(X, path)
        @test model isa X

        m2 = load(path)
        @test m2 isa X

        # TODO test return types here

        S = stoichiometry(model)
        rxns = reactions(model)
        @test n_reactions(model) == length(rxns)
        mets = metabolites(model)
        @test n_metabolites(model) == length(mets)
        gs = genes(model)
        @test n_genes(model) == length(genes)

        @test size(S) == (length(mets), length(rxns))
        @test length(balance(model)) == size(S, 1)
        bs = bounds(model)
        @test bs isa Tuple{Vector{Float64},Vector{Float64}}
        lbs, ubs = bs
        @test length(lbs) == size(S, 2)
        @test length(ubs) == size(S, 2)
        obj = objective(model)
        @test length(obj) == size(S, 2)

        # TODO test that other things work

        # TODO allow (force) folks to add arguments that describe model
        # contents and test their validity right here.

        # TODO test convert
    end
end
