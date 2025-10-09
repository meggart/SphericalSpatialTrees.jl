using SphericalSpatialTrees
using Test

@testset "SphericalSpatialTrees.jl" begin
    include("native.jl")
    include("trees.jl")
    include("test_show_methods.jl")
    include("regulargridtree.jl")
    include("iseatree.jl")
    include("reproject.jl")
end
