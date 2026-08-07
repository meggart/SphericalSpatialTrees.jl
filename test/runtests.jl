using SphericalSpatialTrees
using Test

@testset "SphericalSpatialTrees.jl" begin
    include("nativeisea.jl")
    include("trees.jl")
    include("test_show_methods.jl")
    include("regulargridtree.jl")
    include("utmtree.jl")
    include("iseatree.jl")
    include("equi7tree.jl")
    include("reproject.jl")
    include("lazyprojection.jl")
end
