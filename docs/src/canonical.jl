
# # Implementing new model types
#
# For convenience, `AbstractFBCModels` defines a "canonical" implementation of
# a FBC model: a completely generic data structure that can store exactly the
# complete information that is representable by the `AbstractFBCModel`
# accessors, and nothing else.
#
# The type is not useful for actually constructing the models, but may serve a
# good purpose in several other cases:
#
# - If you need an "intermediate step" when converting complicated FBC models
#   to other types, the canonical model is guaranteed not to lose any
#   information, yet performs relatively well when re-exporting the information
#   via the accessors.
# - If you need to make quick modifications to another model type that does not
#   admit easy mutation (e.g., it is made of immutable `struct`s), you can
#   convert to the canonical model and make the small fixes in there.
# - Here, we use it for describing how to "perform" your own definition of
#   model type, and demonstrate the use of the pre-defined testing framework on
#   it.
#
# The model is available for use as `AbstractFBCModels.CanonicalModel.Model`

# ## Defining the model
#
# For convenience in the later explanation, we list the definition of the
# interface for `AbstractFBCModels.CanonicalModel.Model` here:

# ##LIST src/canonical.jl

# The definition contains several main parts:
#
# - the data structures for the model and all the main parts
# - overloaded accessors that provide generic access for the things in the model
# - overloaded loading and saving functions, together with the declaration of
#   the common model suffixes
# - a conversion function that can extract data using accessors from any other
#   `AbstractFBCModel` and construct the canonical model.
#
# Notably, the default file extension is chosen as very unwieldy so that no one
# ever really exchanges data using this model type.

# ## Testing your model definition
#
# Apart from making sure that the accessors work by usual unit tests, you can
# use 2 testing functions that scrutinize the expected properties of the model
# type both solely as a type, and using an example model file. These allow you
# to discover potential problems, as well as build a self-updating test suite
# for your model that provides long-term sustainability and quality assurance.
#
# ### Running type-level tests
#
# Typically, the test suite would run the following to check if types of
# everything match the expectations.

import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel: Model

A.run_fbcmodel_type_tests(Model);

# ### Making a simple model for value tests
#
# For testing the values, you need to provide an existing file that contains
# the model. Let's create some contents first:

import AbstractFBCModels.CanonicalModel: Reaction, Metabolite, Gene, Coupling

m = Model()
m.metabolites["m1"] = Metabolite(compartment = "inside")
m.metabolites["m2"] = Metabolite(compartment = "outside")
m.genes["g1"] = Gene()
m.genes["g2"] = Gene()
m.reactions["forward"] = Reaction(
    name = "import",
    stoichiometry = Dict("m1" => -1.0, "m2" => 1.0),
    gene_association_dnf = [["g1"], ["g2"]],
    objective_coefficient = 1.0,
)
m.reactions["and_back"] =
    Reaction(name = "export", stoichiometry = Dict("m2" => -1.0, "m1" => 1.0))
m.reactions["exchange1"] = Reaction(
    name = "exchange m1",
    stoichiometry = Dict("m1" => -1.0),
    gene_association_dnf = [[]], # DNF encoding of a reaction that requires no gene products
)
m.reactions["exchange2"] = Reaction(
    name = "exchange m2",
    stoichiometry = Dict("m2" => -1.0),
    gene_association_dnf = [], # DNF encoding of a reaction that never has gene products available
)
m.couplings["total_exchange_limit"] = Coupling(
    lower_bound = 0,
    upper_bound = 10,
    reaction_weights = Dict("exchange$i" => 1.0 for i = 1:2),
)
nothing #hide

show_contains(x, y) = contains(sprint(show, MIME"text/plain"(), x), y) #src
@test show_contains(m, "reactions = ") #src
@test show_contains(m, "metabolites = ") #src
@test show_contains(m, "genes = ") #src
@test show_contains(m.reactions["and_back"], "\"export\"") #src
@test show_contains(m.reactions["forward"], "\"g1\"") #src
@test show_contains(m.metabolites["m1"], "\"inside\"") #src
@test show_contains(m.genes["g1"], "name = nothing") #src
@test show_contains(m.genes["g2"], "name = nothing") #src
@test show_contains(m.couplings["total_exchange_limit"], "name = nothing") #src
@test show_contains(m.couplings["total_exchange_limit"], "upper_bound = 10") #src

# We should immediately find the basic accessors working:
A.stoichiometry(m)

#

A.objective(m)

#

A.coupling(m)

# We can check various side things, such as which reactions would and would not work given all gene products disappear:
products_available = [
    A.reaction_gene_products_available(m, rid, _ -> false) for
    rid in ["forward", "and_back", "exchange1", "exchange2"]
]

@test products_available == [false, nothing, true, false] #src

# We can now also write the model to disk and try to load it with the default
# loading function:
mktempdir() do dir
    path = joinpath(dir, "model.canonical-serialized-fbc")
    A.save(m, path)
    A.load(path)
end

# ### Running the value tests
#
# Given the data, value tests have an opportunity to scrutinize much greater
# amount of properties of the model implementation.
#
# Running the tests requires a model type and an "example" model file:

mktempdir() do dir
    path = joinpath(dir, "model.canonical-serialized-fbc")
    A.save(m, path)
    A.run_fbcmodel_file_tests(Model, path, name = "small model")
end;

# ### Changing identifiers in the model
#
# For convenience, function [`identifier_map`](@ref) can be used to mass-update
# all identifiers and references in a whole canonical model. This is useful
# e.g. when exporting model formats that require certain formatting of the
# identifiers, such as SBML. We demonstrate this by changing all identifiers in
# the above model to have a "type prefix" -- reaction identifiers start with
# `R_`, metabolites start with `M_`, etc.

import AbstractFBCModels.CanonicalModel: identifier_map

m2 = identifier_map(
    m,
    reaction_map = id -> "R_" * id,
    metabolite_map = id -> "M_" * id,
    gene_map = id -> "G_" * id,
    compartment_map = id -> "C_" * id,
    coupling_map = id -> "X_" * id,
);

@test all(startswith("R_"), keys(m2.reactions)) #src
@test all(startswith("M_"), keys(m2.metabolites)) #src
@test all(startswith("G_"), keys(m2.genes)) #src
@test all( #src
    startswith("C_"), #src
    [v.compartment for v in values(m2.metabolites) if !isnothing(v.compartment)], #src
) #src
@test all(startswith("X_"), keys(m2.couplings)) #src

mktempdir() do dir #src
    path = joinpath(dir, "model.canonical-serialized-fbc") #src
    A.save(m2, path) #src
    res = A.run_fbcmodel_file_tests(Model, path, name = "reidentified model") #src
end #src

# We can see that the identifiers have changed:

m2.reactions
@test "R_and_back" in keys(m2.reactions) #src

# ...together with the references:

m2.reactions["R_forward"].stoichiometry

@test "M_m1" in keys(m2.reactions["R_forward"].stoichiometry) #src
