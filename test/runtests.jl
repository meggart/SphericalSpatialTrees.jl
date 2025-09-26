using SphericalSpatialTrees
using Test

@testset "SphericalSpatialTrees.jl" begin
    include("native.jl")
    include("test_show_methods.jl")
end
