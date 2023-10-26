
# # Utility functions
#
# AbstractFBCModels export several utilities for easier definition of the model
# contents; here we explore some of them in individual sections.

import AbstractFBCModels as A

# ## Which accessors should I overload?

# You can query the current set of accessor functions as follows:

A.accessors()

@test length(A.accessors()) == 27 #src

# You do not need to overload all of them (e.g., if you model does not have any
# genes you can completely omit all gene-related functions). The main required
# ones (esp. the reaction- and metabolite-related ones) will throw an error if
# not implemented, and the automated tests will fail with them.

# ## Why are there 2 representations of gene-reaction associations?
#
# Gene-reaction association is, in general, a Boolean function -- you give it a
# list of genes that are available (or, symmetrically, knocked out) and a
# reaction, and it tells you whether the reaction can work or not. This
# representation is captured in accessor `reaction_gene_products_available`.
#
# Some model construction methods, on the other hand, require a relatively
# streamlined view of the gene-reaction associations, roughly describing "ways
# in which several possible enzymes that catalyze the reaction can be assembled
# from individual gene products" -- these are sometimes dubbed "isozymes". This
# view is isomorphic to a Boolean formula in a [disjunctive normal form
# (DNF)](https://en.wikipedia.org/wiki/Disjunctive_normal_form) (put simply, a
# DNF would be "a big OR of ANDed identifiers", such as
# `(a && b && c) || (a && d)`).
# The DNF view of the gene association is accessed via
# `reaction_gene_association_dnf`.
#
# Both representations have pros and cons; the main can be summarized as follows:
# - DNF may not be able to represent some valid boolean formulas without an
#   exponential explosion in size -- in particular, the formula `(a1||b1) &&
#   (a2||b2) && (a3||b3) && ... && (aN||bN)` will explode to roughly `2^N` terms
#   in DNF.
# - There's no simple universal format to store and exchange the descriptions
#   of general boolean functions; we thus do not want to enforce a single one.
#
# Your model is primarily supposed to overload the
# `reaction_gene_products_available` accessor. If you are able to reliably
# convert your Boolean functions to DNF, you should also overload the
# `reaction_gene_association_dnf`.
#
# In some cases when the model directly stores the DNF form, you may also
# utilize the `reaction_gene_products_available_from_dnf` function to easily
# implement the reaction_gene_products_available` by just forwarding the
# arguments.
#
# ## Downloading models for testing
#
# For reproducibility, it is often viable to check downloaded files for
# validity, using e.g. checksums. Since this is a common operation, we provide
# `download_data_file`, which is an appropriate wrapper for
# `Downloads.download`:

mktempdir() do dir
    origin = joinpath(dir, "origin")
    url = "file://$origin"
    dest = joinpath(dir, "model")
    open(origin, "w") do f
        write(f, "hello")
    end
    A.download_data_file(
        url,
        dest,
        "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
    )
    x = A.download_data_file( #src
        url, #src
        dest, #src
        "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824", #src
    ) #src
    @test x == dest #src
    @test read(origin) == read(dest) #src
    open(dest, "w") do f #src
        write(f, "olleh") #src
    end #src
    @test_warn "different" A.download_data_file( #src
        url, #src
        dest, #src
        "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824", #src
    ) #src
end
