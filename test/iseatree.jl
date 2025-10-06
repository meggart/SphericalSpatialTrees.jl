# using Test for unit testing
using Test
import SphericalSpatialTrees as SST
import GeometryOps

@testset "ISEACircleTree construction and properties" begin
	tree = SST.ISEACircleTree(2)
	@test isa(tree, SST.ISEACircleTree)
	@test tree.resolution == 2
	@test length(tree.xr) == 2^2 + 1
	@test length(tree.yr) == 2^2 + 1
	@test SST.gridsize(tree) == (4, 4, 10)
	@test Base.ndims(tree) == 3
	@test SST.nchild(tree) == 10
	@test !SST.isleaf(tree)
	@test SST.nleaf(tree) == 10 * 2^(2*tree.resolution)
	@test isa(SST.extent(tree), GeometryOps.UnitSpherical.SphericalCap)
	@test SST.node_extent(tree) == SST.extent(tree)
	@test SST.rootnode(tree) === tree
	xr, yr = SST.get_xyranges(tree)
	@test xr == tree.xr
	@test yr == tree.yr
end

@testset "ISEACircleTree getchild and index methods" begin
	tree = SST.ISEACircleTree(2)
	# getchild(i) returns a tree node
	child = SST.getchild(tree, 1)
	@test child isa SST.TreeNode
    @test child.grid isa SST.RegularGridTree
	@test child.grid.x == 0:0.25:1.0
    @test child.grid.y == 0:0.25:1.0
    @test child.grid.trans((0.0,1.0)) == [0.0,0.0,1.0]
    @test child.grid.trans((1.0,1.0)) == [-0.8506508083520399, 0.276393202250021, 0.4472135954999579]
    # getchild with vector
	child2 = SST.getchild(tree, [2])
	@test child2 isa SST.TreeNode
	# getchild() returns an iterator of children
	children = collect(SST.getchild(tree))
	@test length(children) == 10
	# index_to_cartesian
	cart = SST.index_to_cartesian(86, tree)
	@test cart == (2,2,6)
	# index_to_unitsphere (should return a point on the sphere)
	pt = SST.index_to_unitsphere(1, tree)
	pt == [-0.1296423624282429, 0.8672401794301382, 0.4807154345826705]
	# index_to_polygon_unitsphere
	poly = SST.index_to_polygon_unitsphere(1, tree)
	@test length(poly) == 5
	poly2 = SST.index_to_polygon_unitsphere(CartesianIndex(1,1,1), tree)
	@test length(poly2) == 5
    @test poly == poly2
    @test poly[1] == child.grid.trans((child.grid.x[1],child.grid.y[1]))
    @test poly[2] == child.grid.trans((child.grid.x[2],child.grid.y[1]))
    @test poly[3] == child.grid.trans((child.grid.x[2],child.grid.y[2]))
    @test poly[4] == child.grid.trans((child.grid.x[1],child.grid.y[2]))
    @test poly[5] == child.grid.trans((child.grid.x[1],child.grid.y[1]))
end

@testset "ISEACircleTree get_gridextent and get_subtree" begin
	tree = SST.ISEACircleTree(2)
	xr, yr = SST.get_xyranges(tree)
    @test xr == 0:0.25:1.0
    @test yr == 0:0.25:1.0
	nr = 1:10
	# get_gridextent returns a SphericalCap
	ext = SST.get_gridextent(tree, 1:2, 1:2, 1:1)
	@test isa(ext, GeometryOps.UnitSpherical.SphericalCap)
	# get_subtree returns a tree node
	subtree = SST.TreeNode(tree, (1:2,1:2,5:5))
	@test subtree isa SST.TreeNode
	@test subtree.index == SST.TreeIndex((1,3),(1,3))
    @test subtree.grid == SST.getchild(tree,5).grid
    subtree2 = SST.TreeNode(tree, (1, 1, [1,2]))
	@test subtree2 === tree
end
