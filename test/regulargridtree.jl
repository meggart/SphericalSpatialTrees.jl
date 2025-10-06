using Test
import SphericalSpatialTrees as SST
using GeometryOps.UnitSpherical: UnitSphereFromGeographic, SphericalCap

@testset "RegularGridTree construction and properties" begin
    x = range(0.0, 1.0, length=5)
    y = range(0.0, 1.0, length=5)
    tree = SST.RegularGridTree(x, y)
    @test isa(tree, SST.RegularGridTree)
    @test length(tree.x) == 5
    @test length(tree.y) == 5
    @test tree.trans == UnitSphereFromGeographic()
    @test SST.gridsize(tree) == (4, 4)
    @test Base.ndims(tree) == 2
    @test SST.nleaf(tree) == 16
    @test SST.get_projection(tree) == tree.trans
    @test SST.get_tag(tree) === nothing
    @test SST.nlevel(tree) == 2
end

@testset "RegularGridTree indexing and geometry" begin
    x = range(0.0, 1.0, length=5)
    y = range(0.0, 1.0, length=5)
    tree = SST.RegularGridTree(x, y)
    # rootnode
    root = SST.rootnode(tree)
    @test isa(root, SST.TreeNode)
    @test root.grid === tree
    @test root.index == SST.TreeIndex((1,5),(1,5))
    # extent
    ext = SST.extent(tree, 1, length(x), 1, length(y))
    @test ext.X == (x[1], x[end])
    @test ext.Y == (y[1], y[end])
    # get_gridextent
    xr = 1:4
    yr = 1:4
    ext2 = SST.get_gridextent(tree, xr, yr)
    @test isa(ext2, SphericalCap)
    # index_to_cartesian
    cart = SST.index_to_cartesian(5, tree)
    @test cart == (1, 2)
    # index_to_unitsphere
    pt = SST.index_to_unitsphere(1, tree)
    @test pt == [0.9999952403603672, 0.0021816546423732855, 0.0021816598343367697]
    # node_to_polygon_unitsphere
    poly = SST.node_to_polygon_unitsphere(root)
    @test length(poly) == 5
    @test poly[1] == root.grid.trans((root.grid.x[1],root.grid.y[1]))
    @test poly[2] == root.grid.trans((root.grid.x[end],root.grid.y[1]))
    @test poly[3] == root.grid.trans((root.grid.x[end],root.grid.y[end]))
    @test poly[4] == root.grid.trans((root.grid.x[1],root.grid.y[end]))
    @test poly[5] == root.grid.trans((root.grid.x[1],root.grid.y[1]))
end

@testset "RegularGridTree TreeNode and children" begin
    x = range(0.0, 1.0, length=5)
    y = range(0.0, 1.0, length=5)
    tree = SST.RegularGridTree(x, y)
    root = SST.rootnode(tree)
    @test SST.nchild(root) == 4
    children = collect(SST.getchild(root))
    @test length(children) == 4
    for c in children
        @test isa(c, SST.TreeNode)
    end
    @test !SST.isleaf(root)
    # Make a leaf node
    leaf = SST.getchild(root, 1)
    while !SST.isleaf(leaf)
        leaf = SST.getchild(leaf, 1)
    end
    @test SST.isleaf(leaf)
    @test isa(SST.node_extent(leaf), SphericalCap)
    @test SST.linind(leaf) == 1
end

