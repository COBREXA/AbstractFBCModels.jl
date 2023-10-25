
using Test

import AbstractFBCModels as A

@testset "AbstractFBCModels tests" begin
    @testset "Canonical model documentation" begin
        include("../docs/src/canonical.jl")
    end
    @testset "Utilities" begin
        include("../docs/src/canonical.jl")
    end
end
