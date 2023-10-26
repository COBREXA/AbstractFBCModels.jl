
# # Custom models
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
#   information, yet perform relatively well when re-exporting the information
#   via the accessors.
# - If you need to make quick modifications to another model type that does not
#   admin easy mutation (e.g., it is made of immutable `struct`s), you can
#   convert to the canonical model and make the small fixes in there.
# - Here, we use it for describing how to "perform" your own definition of
#   model type, and demonstrate the use of the pre-defined testing framework on
#   it.
#
# The model is available for use as `AbstractFBCModels.CanonicalModel.Model`

# ## Defining the model
#
# For convenience in the later explanation, we list the whole definition of the
# module here:

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

import AbstractFBCModels.CanonicalModel: Reaction, Metabolite, Gene

m = Model()
m.metabolites["m1"] = Metabolite(compartment = "inside")
m.metabolites["m2"] = Metabolite(compartment = "outside")
m.genes["g1"] = Gene()
m.genes["g2"] = Gene()
m.reactions["forward"] = Reaction(
    name = "import",
    stoichiometry = Dict("m1" => -1.0, "m2" => 1.0),
    gene_association_dnf = [["g1"], ["g2"]],
)
m.reactions["and_back"] =
    Reaction(name = "export", stoichiometry = Dict("m2" => -1.0, "m1" => 1.0))
nothing #hide

# We should immediately find the basic accessors working:
A.stoichiometry(m)

# We can check various side things, such as which reactions would and would not work given all gene products disappear:
products_available = [
    A.reaction_gene_products_available(m, rid, _ -> false) for
    rid in ["forward", "and_back"]
]

@test products_available == [false, nothing] #src

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
