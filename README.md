
# AbstractFBCModels.jl -- Common interface for flux-balanced constrained models

| Build status | Documentation |
|:---:|:---:|
| [![CI](https://github.com/COBREXA/AbstractFBCModels.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/COBREXA/AbstractFBCModels.jl/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/COBREXA/AbstractFBCModels.jl/branch/master/graph/badge.svg?token=A2ui7exGIH)](https://codecov.io/gh/COBREXA/AbstractFBCModels.jl) | [![stable documentation](https://img.shields.io/badge/docs-stable-blue)](https://cobrexa.github.io/AbstractFBCModels.jl/stable) [![dev documentation](https://img.shields.io/badge/docs-dev-cyan)](https://cobrexa.github.io/AbstractFBCModels.jl/dev) |

Package `AbstractFBCModels.jl` defines a set of accessor functions that provide
a generic API for working with the contents of the constraint-based metabolic
(or "Flux-Balanced Constrained" for FBC) models. Packages that work with model
formats may implement the API to transparently expose the contents of the
models in given format to many other packages.

The primary purpose of this is to provide the model loading functionality for
[COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl) and
[FBCModelTests.jl](https://github.com/LCSB-BioCore/FBCModelTests.jl), but is
otherwise completely generic and can be used independently of these packages.

The package is currently quite new, maintained and open for extensions; feel
free to discuss changes and ideas via issues and pull requests.

#### Acknowledgements

`AbstractFBCModels.jl` was developed at the Luxembourg Centre for Systems
Biomedicine of the University of Luxembourg
([uni.lu/lcsb](https://www.uni.lu/lcsb))
and at Institute for Quantitative and Theoretical Biology at Heinrich Heine
University Düsseldorf ([qtb.hhu.de](https://www.qtb.hhu.de/en/)).
The development was supported by European Union's Horizon 2020 Programme under
PerMedCoE project ([permedcoe.eu](https://www.permedcoe.eu/)),
agreement no. 951773.

<img src="docs/src/assets/unilu.svg" alt="Uni.lu logo" height="64px">   <img src="docs/src/assets/lcsb.svg" alt="LCSB logo" height="64px">   <img src="docs/src/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto">   <img src="docs/src/assets/permedcoe.svg" alt="PerMedCoE logo" height="64px">
