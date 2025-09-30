using SphericalSpatialTrees
using Test

@testset "SphericalSpatialTrees.jl" begin
    include("native.jl")
    include("trees.jl")
    include("test_show_methods.jl")
end
