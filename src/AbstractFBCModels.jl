"""
Package `AbstractFBCModels.jl` defines a very general interface that allows
access to data stored in flux-balanced constrained (FBC) models, which is a
common way to represent metabolic models of organisms of various scales.

The interface enables 2 main goals:
- You can load and save data from/to various types of FBC models, including the
  MatLab-based format (as used e.g. by
  [CobraToolbox](https://github.com/opencobra/cobratoolbox)), JSON-based format
  (adopted by [CobraPy](https://github.com/opencobra/cobrapy/)),
  [SBML](https://github.com/sbmlteam/libsbml), and others.
- You can freely convert the model data among the formats using standard Julia
  `convert()`.

The package provides an additional "canonical model" format (in submodule
`CanonicalModel`) that implements the bare minimum of features required to
store all data representable via the general interface. You can use it as a
base for implementing more complex model formats, or as a safe middle-point for
model data conversion. See the examples in the documentation for details.

This package is lightweight and implements no other specific functionality. To
load data from actual model formats, you will also need other packages that
implement the functionality, such as `SBMLFBCModels.jl`.
"""
module AbstractFBCModels

using DocStringExtensions

include("types.jl")
include("accessors.jl")
include("io.jl")
include("utils.jl")
include("canonical.jl")

end # module AbstractFBCModel
