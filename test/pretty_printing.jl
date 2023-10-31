@testset "Canonical types printing" begin
    # create and test IO for dummy gene
    g = A.CanonicalModel.Gene()
    g.name = "gene1"
    g.notes = Dict("notes" => ["blah", "blah"])
    g.annotations = Dict("sboterm" => ["sbo"], "ncbigene" => ["ads", "asds"])
    @test all(contains.(sprint(show, MIME("text/plain"), g), ["gene1", "blah", "asds"]))

    # create and test IO for dummy reaction
    r1 = A.CanonicalModel.Reaction()
    r1.name = "r1"
    r1.stoichiometry = Dict("m1" => -1.0, "m2" => 1.0)
    r1.lower_bound = -100.0
    r1.upper_bound = 100.0
    r1.gene_association_dnf = [["g1", "g2"], ["g3"]]
    r1.notes = Dict("notes" => ["blah", "blah"])
    r1.annotations = Dict("sboterm" => ["sbo"], "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), r1),
            ["r1", "100.0", "↔", "(g1 and g2)", "(g3)", "blah", "biocyc"],
        ),
    )

    r1.stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => -1.0, "m4" => -1.0, "m5" => -1.0, "m6" => -1.0,)
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), r1),
            ["-100.0", "∅", "..."],
        ),
    )

    r1.lower_bound = -100.0
    r1.upper_bound = 0.0
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), r1),
            ["←"],
        ),
    )

    r1.lower_bound = 0.0
    r1.upper_bound = 100.0
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), r1),
            ["→"],
        ),
    )

    r1.lower_bound = 0.0
    r1.upper_bound = 0.0
    @test all(
        contains.(
            sprint(show, MIME("text/plain"), r1),
            ["→|←"],
        ),
    )

    # create and test IO for dummy metabolite
    m1 = A.CanonicalModel.Metabolite()
    m1.name = "met1"
    m1.formula = Dict("C" => 6, "H" => 12, "O" => 6)
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes" => ["blah", "blah"])
    m1.annotations = Dict("sboterm" => ["sbo"], "kegg.compound" => ["ads", "asds"])

    @test all(
        contains.(
            sprint(show, MIME("text/plain"), m1),
            ["met1", "C6H12O6", "blah", "asds"],
        ),
    )
end