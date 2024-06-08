
@testset "behavior of accessors on abstract type" begin

    import AbstractFBCModels as A

    struct NotAModel <: A.AbstractFBCModel end

    m = NotAModel()

    @test_throws ErrorException A.reactions(m)
    @test_throws ErrorException A.metabolites(m)
    @test_throws ErrorException A.genes(m)
    @test_throws ErrorException A.n_reactions(m)
    @test_throws ErrorException A.n_metabolites(m)
    @test_throws ErrorException A.n_genes(m)
    @test_throws ErrorException A.stoichiometry(m)
    @test_throws ErrorException A.balance(m)
    @test_throws ErrorException A.bounds(m)
    @test_throws ErrorException A.objective(m)

    @test isnothing(A.reaction_gene_products_available(m, "", _ -> True))
    @test isnothing(A.reaction_gene_association_dnf(m, ""))
    @test_throws ErrorException A.reaction_stoichiometry(m, "")
    @test isnothing(A.metabolite_formula(m, ""))
    @test isnothing(A.metabolite_charge(m, ""))
    @test isnothing(A.metabolite_compartment(m, ""))
    @test isempty(A.couplings(m))
    @test A.a_couplings(m) == 0
    @test_throws ErrorException A.coupling(m)
    @test all(isempty, A.coupling_bounds(m))
    @test_throws ErrorException A.coupling_weights(m, "")

    @test isempty(A.reaction_annotations(m, ""))
    @test isempty(A.metabolite_annotations(m, ""))
    @test isempty(A.gene_annotations(m, ""))
    @test isempty(A.coupling_annotations(m, ""))
    @test isempty(A.reaction_notes(m, ""))
    @test isempty(A.metabolite_notes(m, ""))
    @test isempty(A.gene_notes(m, ""))
    @test isempty(A.coupling_notes(m, ""))
    @test isnothing(A.reaction_name(m, ""))
    @test isnothing(A.metabolite_name(m, ""))
    @test isnothing(A.gene_name(m, ""))
    @test isnothing(A.coupling_name(m, ""))

    @test_throws ErrorException A.load(NotAModel, ".")
    @test_throws ErrorException A.save(m, ".")
    @test_throws ErrorException A.filename_extensions(NotAModel)
    @test_throws ErrorException show(stdout, MIME"text/plain"(), m)

    A.n_reactions(::NotAModel) = 123

    @test size(A.coupling(m) == (0, 123))
end
