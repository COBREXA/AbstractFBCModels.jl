
import Test
import Test: @testset, @test

"""
$(TYPEDSIGNATURES)

Internal helper macro for rendering easily interpretable test failures.
"""
macro atest(ex, nice, actually = nothing)
    result = quote
        try
            $(Test.Returned)($(esc(ex)), $(esc(actually)), $(QuoteNode(__source__)))
        catch _e
            _e isa InterruptException && rethrow()
            $(Test.Threw)(_e, Base.current_exceptions(), $(QuoteNode(__source__)))
        end
    end
    Base.remove_linenums!(result)
    Expr(:block, __source__, Expr(:call, Test.do_test, result, esc(nice)))
end

"""
$(TYPEDSIGNATURES)

Test if the given model type loads properly from a file.

The function uses the testing infrastructure from `Test` to report problems --
it is supposed to be a part of larger test-sets, preferably in all model
implementation packages.
"""
function run_fbcmodel_file_tests(
    ::Type{X},
    path::String;
    name::String = path,
    test_save::Bool = true,
) where {X<:AbstractFBCModel}
    @testset "Model `$name' of type $X" begin

        model = load(X, path)
        @test model isa X

        m2 = load(path)
        @test m2 isa X

        S = stoichiometry(model)
        C = coupling(model)
        rxns = reactions(model)
        mets = metabolites(model)
        gens = genes(model)
        cpls = couplings(model)

        @test n_reactions(model) == length(rxns)
        @test n_metabolites(model) == length(mets)
        @test n_genes(model) == length(gens)
        @test n_couplings(model) == length(cpls)

        # test sizing
        @test size(S) == (length(mets), length(rxns))
        @test size(C) == (length(cpls), length(rxns))
        @test length(balance(model)) == size(S, 1)

        bs = bounds(model)
        @test bs isa Tuple{Vector{Float64},Vector{Float64}}
        lbs, ubs = bs
        @test length(lbs) == size(S, 2)
        @test length(ubs) == size(S, 2)

        cbs = coupling_bounds(model)
        @test cbs isa Tuple{Vector{Float64},Vector{Float64}}
        clbs, cubs = cbs
        @test length(clbs) == size(C, 1)
        @test length(cubs) == size(C, 1)

        obj = objective(model)
        @test length(obj) == size(S, 2)

        let ms = Set(mets), mi = Dict(mets .=> 1:length(mets))
            for (ridx, rid) in enumerate(reactions(model))
                for (met, stoi) in reaction_stoichiometry(model, rid)
                    # test if reaction stoichiometries refer to metabolites and the result is the same as with the matrix
                    @atest met in ms "metabolite `$met' in reaction_stoichiometry() of `$rid' is in metabolites()"
                    @atest S[mi[met], ridx] == stoi "reaction_stoichiometry() of reaction `$rid' matches the column in stoichiometry() matrix"
                end
                # TODO also test the other direction
            end
        end

        let rs = Set(rxns), ri = Dict(rxns .=> 1:length(rxns))
            for (cidx, cid) in enumerate(cpls)
                for (rxn, w) in coupling_weights(model, cid)
                    # test if coupling weights are the same as with the matrix
                    @atest rxn in rs "reaction `$rxn' in coupling_weights() of `$cid' is in reactions()"
                    @atest C[cidx, ri[rxn]] == w "coupling_weights() of coupling `$cid' matches the row in coupling() matrix"
                end
                # TODO also test the other direction
            end
        end

        let gs = Set(gens)
            for rid in rxns
                # check if all genes reported in DNFs are registered in genes
                ga = reaction_gene_association_dnf(model, rid)
                isnothing(ga) && continue
                gas = vcat(ga...)
                @atest all(g -> g in gs, gas) "all gene products referenced by DNF of `$rid' are in genes()"
                # if we have a DNF for some reaction, this must be able to interpret it
                @atest !isnothing(reaction_gene_products_available(model, rid, _ -> true)) "reaction_gene_products_available() returns a Boolean for reaction `$rid' with known DNF"
                # ...and the interpretation must be the same.
                if isempty(ga)
                    @atest !reaction_gene_products_available(
                        model,
                        rid,
                        g -> error("should not query"),
                    ) "reaction_gene_products_available(#= unsatisfiable DNF of `$rid' =#) == false"
                elseif [] in ga
                    @atest reaction_gene_products_available(model, rid, g -> false) "reaction_gene_products_available(#= tautological DNF of `$rid' =#) == true"
                else
                    # check if the DNF descriptions match what the boolean evaluators think
                    sgas = Set(gas)
                    @atest reaction_gene_products_available(model, rid, g -> g in sgas) "reaction_gene_products_available(#= all DNF terms for `$rid' are true =#) == true"
                    @atest !reaction_gene_products_available(model, rid, g -> !(g in sgas)) "reaction_gene_products_available(#= all DNF terms for `$rid' are false =#) == false"
                end
            end
        end

        # test saving
        m3 = test_save ? mktempdir() do d
            p = joinpath(d, "test-model")
            save(model, p)
            load(X, p)
        end : m2

        # test conversion through the canonical model (should preserve everything)
        m4 = convert(X, convert(CanonicalModel.Model, model))

        # all the variants loaded from other sources should be equivalent
        for m in [m2, m3, m4]
            @test issetequal(rxns, reactions(m))
            @test issetequal(mets, metabolites(m))
            @test issetequal(gens, genes(m))
            @test issetequal(cpls, couplings(m))

            @test Dict(rxns .=> collect(obj)) ==
                  Dict(reactions(m) .=> collect(objective(m)))

            for rxn in rxns
                @test issetequal(reaction_subsystem(model, rxn), reaction_subsystem(m, rxn))
            end

            # TODO eventually add stricter equality tests
        end

        # TODO allow users to add arguments that describe model contents (like
        # expect_n_reactions=123) and test their validity right here.
    end
end

export run_fbcmodel_file_tests

"""
$(TYPEDSIGNATURES)

Test if the given model type works right.

The function uses the testing infrastructure from `Test` to report problems --
it is supposed to be a part of larger test-sets, preferably in all model
implementation packages.
"""
function run_fbcmodel_type_tests(::Type{X}) where {X<:AbstractFBCModel}
    @testset "Model type $X properties" begin
        rt(f, ret, args...) = @testset "$f should return $ret" begin
            @atest all(t -> t <: ret, unique(Base.return_types(f, args))) "all return types of $f$args are subtypes of $ret"
        end

        rt(reactions, Vector{String}, X)
        rt(n_reactions, Int, X)
        rt(metabolites, Vector{String}, X)
        rt(n_metabolites, Int, X)
        rt(genes, Vector{String}, X)
        rt(n_genes, Int, X)
        rt(couplings, Vector{String}, X)
        rt(n_couplings, Int, X)
        rt(stoichiometry, SparseMat, X)
        rt(coupling, SparseMat, X)
        rt(bounds, Tuple{Vector{Float64},Vector{Float64}}, X)
        rt(coupling_bounds, Tuple{Vector{Float64},Vector{Float64}}, X)
        rt(objective, SparseVec, X)
        rt(balance, SparseVec, X)

        rt(reaction_stoichiometry, Dict{String,Float64}, X, String)
        rt(reaction_gene_products_available, Maybe{Bool}, X, String, Function)
        rt(reaction_gene_association_dnf, Maybe{GeneAssociationDNF}, X, String)
        rt(reaction_name, Maybe{String}, X, String)
        rt(reaction_annotations, Annotations, X, String)
        rt(reaction_notes, Notes, X, String)

        rt(metabolite_formula, Maybe{MetaboliteFormula}, X, String)
        rt(metabolite_charge, Maybe{Int}, X, String)
        rt(metabolite_compartment, Maybe{String}, X, String)
        rt(metabolite_name, Maybe{String}, X, String)
        rt(metabolite_annotations, Annotations, X, String)
        rt(metabolite_notes, Notes, X, String)

        rt(gene_name, Maybe{String}, X, String)
        rt(gene_annotations, Annotations, X, String)
        rt(gene_notes, Notes, X, String)

        rt(coupling_weights, Dict{String,Float64}, X, String)
        rt(coupling_name, Maybe{String}, X, String)
        rt(coupling_annotations, Annotations, X, String)
        rt(coupling_notes, Notes, X, String)
    end
end

export run_fbcmodel_type_tests
