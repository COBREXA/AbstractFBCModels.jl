
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

# The definition contains 3 main parts:
#
# - the data structures for the model and all the main parts
# - overloaded accessors that provide generic access for the things in the model
# - overloaded loading and saving functions, together with the declaration of
#   the common model suffixes
# - a conversion function that can extract data using accessors from any other
#   `AbstractFBCModel` and construct the canonical model.

# ## Testing your model definition
#
# Apart from making sure that the accessors work by usual unit tests, you can
# use 2 testing functions that scrutinize the expected properties of the model
# type both solely as a type, and using an example model file. These allow you
# to discover potential problems, as well as build a self-updating test suite
# for your model that provides long-term sustainability and quality assurance.
#
# Typically, in your test suite, you would run these functions:

import AbstractFBCModels as A

A.run_fbcmodel_type_tests(A.CanonicalModel.Model);
