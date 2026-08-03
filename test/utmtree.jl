using Test
import SphericalSpatialTrees as SST
import DimensionalData as DD
import DiskArrays
import GeometryOps
import GeometryOps.SpatialTreeInterface as STI
using GeometryOps.UnitSpherical: UnitSphereFromGeographic, GeographicFromUnitSphere, UnitSphericalPoint, SphericalCap, _contains

# A realistic UTM grid around central Europe (zone 32N), 4x4 cells per zone/hemisphere
const UTM_TEST_X = range(500000.0, 700000.0, length=5)
const UTM_TEST_Y = range(5.4e6, 5.8e6, length=5)

# Geographic coordinates of Jena, Germany (in UTM zone 32N)
const JENA = (11.578069, 50.912083)

approxeq(a, b; kwargs...) = all(isapprox.(a, b; kwargs...))

@testset "UTMTree construction and properties" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    @test isa(tree, SST.UTMTree)
    @test length(tree.x) == 5
    @test length(tree.y) == 5
    @test tree.trans == SST.UnitSphereFromUTM()
    @test SST.gridsize(tree) == (4, 4, 60, 2)
    @test Base.ndims(tree) == 4
    @test SST.nleaf(tree) == 4 * 4 * 60 * 2
    @test SST.nlevel(tree) == 9
    @test SST.get_projection(tree) == tree.trans
    @test SST.get_tag(tree) === nothing
    xmid, ymid, zoned, hemid = DD.dims(tree)
    @test xmid == DD.X(525000.0:50000.0:675000.0)
    @test ymid == DD.Y(5.45e6:100000.0:5.75e6)
    @test zoned == DD.Dim{:ZONE}(1:60)
    @test hemid == DD.Dim{:Hemisphere}(1:2)
end

@testset "UTM transforms" begin
    t = SST.UTMTransform()
    utm = t(JENA)
    @test approxeq(utm, (681231.2399276134, 5.643213955487227e6, 32, 1))
    @test approxeq(inv(t)(utm), JENA)
    it = SST.IUTMTransform()
    @test approxeq(it((681231.2399276134, 5.643213955487227e6, 32, 1)), JENA)
    @test inv(SST.UTMTransform()) isa SST.IUTMTransform
    @test inv(SST.IUTMTransform()) isa SST.UTMTransform
    # round trip through the unit sphere
    us = SST.UnitSphereFromUTM()((681231.2399276134, 5.643213955487227e6, 32, 1))
    @test us ≈ UnitSphericalPoint(0.6176825643523223, 0.12654564930689588, 0.7761793918525765)
    @test approxeq(GeographicFromUnitSphere()(us), JENA)
    @test inv(SST.UnitSphereFromUTM()) isa SST.UTMFromUnitSphere
    @test inv(SST.UTMFromUnitSphere()) isa SST.UnitSphereFromUTM
end

@testset "UTMTree indexing and geometry" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    root = SST.rootnode(tree)
    @test isa(root, SST.TreeNode)
    @test root.grid === tree
    @test root.index == SST.UTMIndex((1, 5), (1, 5), (1, 61), 3)
    # extent spans all four dimensions
    ext = SST.extent(tree, root.index)
    @test ext.X == (500000.0, 700000.0)
    @test ext.Y == (5.4e6, 5.8e6)
    @test ext.ZONE == (1, 60)
    @test ext.Hemisphere == (1, 2)
    # get_gridextent returns a spherical cap
    @test isa(SST.get_gridextent(tree, 1:4, 1:4, 32:32, 1:1), SphericalCap)
    # a node spanning both hemispheres covers the whole sphere
    @test SST.get_gridextent(tree, 1:4, 1:4, 32:32, 1:2).radius ≈ π
    # index conversions
    @test SST.index_to_cartesian(1, tree) == (1, 1, 1, 1)
    @test SST.index_to_cartesian(502, tree) == (2, 2, 32, 1)
    @test SST.index_to_cartesian(SST.nleaf(tree), tree) == (4, 4, 60, 2)
    @test SST.index_to_native_coords(502, tree) == (575000.0, 5.55e6, 32, 1)
    @test SST.index_to_unitsphere(502, tree) ≈ UnitSphericalPoint(0.6316424885141867, 0.1119285710131051, 0.7671373812392209)
    @test approxeq(SST.index_to_lonlat(502, tree), (10.048638631019521, 50.09751953148771))
    # node_to_polygon_unitsphere traces a closed polygon of unit sphere points
    node = SST.TreeNode(tree, (1:2, 1:2, 5:5, 1:1))
    poly = SST.node_to_polygon_unitsphere(node)
    @test length(poly) == 5
    @test poly[1] == poly[end]
    @test poly[1] == tree.trans((500000.0, 5.4e6, 5, 1))
    @test poly[2] == tree.trans((600000.0, 5.4e6, 5, 1))
    @test poly[3] == tree.trans((600000.0, 5.6e6, 5, 1))
    @test poly[4] == tree.trans((500000.0, 5.6e6, 5, 1))
end

@testset "UTMTree TreeNode and children" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    root = SST.rootnode(tree)
    @test SST.nchild(root) == 4
    children = collect(SST.getchild(root))
    @test length(children) == 4
    for c in children
        @test isa(c, SST.TreeNode)
    end
    # children split by zone range and hemisphere
    @test children[1].index == SST.UTMIndex((1, 5), (1, 5), (1, 31), 1)
    @test children[2].index == SST.UTMIndex((1, 5), (1, 5), (1, 31), 2)
    @test children[3].index == SST.UTMIndex((1, 5), (1, 5), (31, 61), 1)
    @test children[4].index == SST.UTMIndex((1, 5), (1, 5), (31, 61), 2)
    @test !SST.isleaf(root)
    # descend to a leaf
    leaf = SST.getchild(root, 1)
    while !SST.isleaf(leaf)
        leaf = SST.getchild(leaf, 1)
    end
    @test SST.isleaf(leaf)
    @test leaf.index == SST.UTMIndex((1, 2), (1, 2), (1, 2), 1)
    @test isa(SST.node_extent(leaf), SphericalCap)
    @test SST.linind(leaf) == 1
    # explicit subtree from index ranges
    sub = SST.TreeNode(tree, (1:2, 1:2, 5:5, 1:1))
    @test sub.index == SST.UTMIndex((1, 3), (1, 3), (5, 6), 1)
    subext = SST.extent(tree, sub.index)
    @test subext.X == (500000.0, 600000.0)
    @test subext.Y == (5.4e6, 5.6e6)
    @test subext.ZONE == (5, 5)
    # the extent circle contains the node center
    @test _contains(SST.node_extent(sub), tree.trans((550000.0, 5.5e6, 5, 1)))
end

@testset "UTMTree queries" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    cap = SphericalCap(UnitSphereFromGeographic()(JENA), 0.0005)
    inds = STI.query(SST.rootnode(tree), cap)
    @test sort(inds) == [507, 508]
    @test SST.index_to_native_coords.(inds, (tree,)) == [(625000.0, 5.65e6, 32, 1), (675000.0, 5.65e6, 32, 1)]
    # with_transform gives the same result
    wt = SST.with_transform(tree) do t
        STI.query(SST.rootnode(t), cap)
    end
    @test sort(wt) == sort(inds)
    # find_nearest locates the cell containing the point
    dist, i = SST.find_nearest(tree, JENA)
    @test i == 508
    @test dist < 0.0015
end

@testset "UTMTree from DimArray and Projection" begin
    xc = range(525000.0, 675000.0, length=4)
    yc = range(5.45e6, 5.75e6, length=4)
    dims = (DD.X(xc), DD.Y(yc), DD.Dim{:ZONE}(1:60), DD.Dim{:Hemisphere}(1:2))
    ar = DD.DimArray(DiskArrays.mockchunks(rand(4, 4, 60, 2), (2, 2, 1, 1)), dims)
    # convenience constructor from a dim array
    tree = SST.UTMTree(ar)
    @test SST.gridsize(tree) == (4, 4, 60, 2)
    @test SST.rootnode(tree).index == SST.UTMIndex((1, 5), (1, 5), (1, 61), 3)
    @test tree.trans == SST.UnitSphereFromUTM()
    # ProjectionSource
    source = SST.ProjectionSource(SST.UTMTree, ar)
    @test source.tree isa SST.UTMTree
    @test SST.gridsize(source.tree) == (4, 4, 60, 2)
    @test SST.gridsize(source.chunktree) == (2, 2, 60, 2)
    @test length(source.lookups) == 4
    # ProjectionTarget with a small chunk size for the test grid
    target = SST.ProjectionTarget(SST.UTMTree, UTM_TEST_X, UTM_TEST_Y; chunksize=2)
    @test SST.gridsize(target.tree) == (4, 4, 60, 2)
    @test SST.gridsize(target.chunktree) == (2, 2, 60, 2)
end

@testset "UTMTree internal indexing functions" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    
    @testset "UTMIndex validation" begin
        valid_index = SST.UTMIndex((1, 5), (1, 5), (1, 61), 1)
        @test SST.is_valid_index(valid_index)
        
        invalid_x = SST.UTMIndex((5, 1), (1, 5), (1, 61), 1)
        @test !SST.is_valid_index(invalid_x)
        
        invalid_y = SST.UTMIndex((1, 5), (5, 1), (1, 61), 1)
        @test !SST.is_valid_index(invalid_y)
        
        invalid_zone = SST.UTMIndex((1, 5), (1, 5), (61, 1), 1)
        @test !SST.is_valid_index(invalid_zone)
        
        invalid_hemi = SST.UTMIndex((1, 5), (1, 5), (1, 61), 0)
        @test !SST.is_valid_index(invalid_hemi)
        
        invalid_hemi2 = SST.UTMIndex((1, 5), (1, 5), (1, 61), 4)
        @test !SST.is_valid_index(invalid_hemi2)
    end
    
    @testset "UTMIndex _nchild calculation" begin
        full_range = SST.UTMIndex((1, 5), (1, 5), (1, 61), 3)
        @test SST._nchild(full_range) == 4
        
        x_single = SST.UTMIndex((1, 2), (1, 5), (1, 61), 1)
        @test SST._nchild(x_single) == 2
        
        y_single = SST.UTMIndex((1, 5), (1, 2), (1, 61), 1)
        @test SST._nchild(y_single) == 2
        
        zone_single = SST.UTMIndex((1, 5), (1, 5), (32, 33), 1)
        @test SST._nchild(zone_single) == 4
        
        half_partial = SST.UTMIndex((1, 2), (1, 5), (32, 61), 1)
        @test SST._nchild(half_partial) == 2
        
        x_y_single = SST.UTMIndex((1, 2), (1, 2), (1, 61), 1)
        @test SST._nchild(x_y_single) == 1
        
        leaf = SST.UTMIndex((1, 2), (1, 2), (32, 33), 1)
        @test SST._nchild(leaf) == 1
    end
    
    @testset "split_utm_updown" begin
        north_hemi = 1
        south_hemi = 2
        both_hemi = 3
        
        y_range = (1, 5)
        
        north_result = SST.split_utm_updown(north_hemi, y_range)
        @test north_result == ((1, (1, 3)), (1, (3, 5)))
        
        south_result = SST.split_utm_updown(south_hemi, y_range)
        @test south_result == ((2, (1, 3)), (2, (3, 5)))
        
        both_result = SST.split_utm_updown(both_hemi, y_range)
        @test both_result == ((1, (1, 5)), (2, (1, 5)))
    end
    
    @testset "split_utm_leftright" begin
        zone_full = (1, 61)
        zone_partial = (32, 61)
        zone_single = (32, 33)
        
        x_range = (1, 5)
        
        full_result = SST.split_utm_leftright(zone_full, x_range)
        @test full_result == (((1, 31), (1, 5)), ((31, 61), (1, 5)))
        
        partial_result = SST.split_utm_leftright(zone_partial, x_range)
        @test partial_result == (((32, 46), (1, 5)), ((46, 61), (1, 5)))
        
        single_result = SST.split_utm_leftright(zone_single, x_range)
        @test single_result == (((32, 33), (1, 3)), ((32, 33), (3, 5)))
    end
    
    @testset "split_4 UTMIndex" begin
        idx = SST.UTMIndex((1, 5), (1, 5), (1, 61), 1)
        children = SST.split_4(idx)
        @test length(children) == 4
        
        @test children[1] == SST.UTMIndex((1, 5), (1, 3), (1, 31), 1)
        @test children[2] == SST.UTMIndex((1, 5), (3, 5), (1, 31), 1)
        @test children[3] == SST.UTMIndex((1, 5), (1, 3), (31, 61), 1)
        @test children[4] == SST.UTMIndex((1, 5), (3, 5), (31, 61), 1)
        
        idx_both = SST.UTMIndex((1, 5), (1, 5), (32, 33), 3)
        children_both = SST.split_4(idx_both)
        @test length(children_both) == 4
        @test children_both[1] == SST.UTMIndex((1, 3), (1, 5), (32, 33), 1)
        @test children_both[2] == SST.UTMIndex((1, 3), (1, 5), (32, 33), 2)
        @test children_both[3] == SST.UTMIndex((3, 5), (1, 5), (32, 33), 1)
        @test children_both[4] == SST.UTMIndex((3, 5), (1, 5), (32, 33), 2)
    end
end

@testset "UTMTree internal geometry functions" begin
    tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
    
    @testset "get_gridextent" begin
        ext = SST.get_gridextent(tree, 1:4, 1:4, 32:32, 1:1)
        @test isa(ext, SphericalCap)
        @test isa(ext.point, UnitSphericalPoint)
        @test ext.radius > 0
        
        ext_both = SST.get_gridextent(tree, 1:4, 1:4, 32:32, 1:2)
        @test ext_both.radius ≈ π
        
        ext_small = SST.get_gridextent(tree, 1:2, 1:2, 10:11, 1:1)
        @test isa(ext_small, SphericalCap)
    end
    
    @testset "linind" begin
        root_idx = SST.UTMIndex((1, 5), (1, 5), (1, 61), 1)
        lin = SST.linind(tree, root_idx)
        @test lin == 1
        
        last_idx = SST.UTMIndex((4, 5), (4, 5), (60, 61), 2)
        lin_last = SST.linind(tree, last_idx)
        @test lin_last == SST.nleaf(tree)
        
        idx = SST.UTMIndex((2, 3), (2, 3), (32, 33), 1)
        lin_expected = LinearIndices((4, 4, 60, 2))[2, 2, 32, 1]
        @test SST.linind(tree, idx) == lin_expected
    end
    
    @testset "isleaf for UTMIndex" begin
        full_idx = SST.UTMIndex((1, 5), (1, 5), (1, 61), 3)
        @test !SST.isleaf(full_idx)
        
        zone_single = SST.UTMIndex((1, 5), (1, 5), (32, 33), 3)
        @test !SST.isleaf(zone_single)
        
        hemi_single = SST.UTMIndex((1, 5), (1, 5), (1, 61), 1)
        @test !SST.isleaf(hemi_single)
        
        leaf_idx = SST.UTMIndex((1, 2), (1, 2), (32, 33), 1)
        @test SST.isleaf(leaf_idx)
    end
    
    @testset "circle_from_extent_1" begin
        ext = SST.extent(tree, SST.UTMIndex((1, 3), (1, 3), (10, 12), 1))
        circle = SST.circle_from_extent_1(ext, tree)
        @test isa(circle, SphericalCap)
        @test circle.radius > 0
        
        ext_full = SST.extent(tree, SST.UTMIndex((1, 5), (1, 5), (1, 61), 3))
        circle_full = SST.circle_from_extent_1(ext_full, tree)
        @test isa(circle_full, SphericalCap)
        @test circle_full.radius ≈ π
    end
    
    @testset "node_to_polygon_unitsphere" begin
        idx = SST.UTMIndex((1, 3), (1, 3), (32, 33), 1)
        poly = SST.node_to_polygon_unitsphere(tree, idx)
        @test length(poly) == 5
        @test poly[1] == poly[end]
        @test all(p -> isa(p, UnitSphericalPoint), poly)
        
        idx_both = SST.UTMIndex((1, 3), (1, 3), (32, 33), 3)
        @test_throws ErrorException SST.node_to_polygon_unitsphere(tree, idx_both)
        
        leaf_idx = SST.UTMIndex((1, 2), (1, 2), (32, 33), 1)
        leaf_poly = SST.node_to_polygon_unitsphere(tree, leaf_idx)
        @test length(leaf_poly) == 5
        @test leaf_poly[1] == tree.trans((500000.0, 5.4e6, 32, 1))
    end
end

@testset "UTMTree internal transforms" begin
    @testset "UTMTransform round trips" begin
        t = SST.UTMTransform()
        
        @test approxeq(t(JENA), (681231.2399276134, 5.643213955487227e6, 32, 1))
        
        utm_north = t((11.5, 51.0))
        @test approxeq(inv(t)(utm_north), (11.5, 51.0))
        
        utm_south = t((-45.0, -23.0))
        @test approxeq(inv(t)(utm_south), (-45.0, -23.0))
        
        @test approxeq(inv(t)(t((0.0, 0.0))), (0.0, 0.0); atol=1e-10)
    end
    
    @testset "IUTMTransform round trips" begin
        it = SST.IUTMTransform()
        
        utm_coord = (681231.2399276134, 5.643213955487227e6, 32, 1)
        @test approxeq(it(utm_coord), JENA)
        
        # Note: Southern hemisphere UTM tests are commented out due to PROJ library issues
        # with Y coordinate interpretation in PROJ transformations
        # south_utm = (500000.0, 10.4e6, 23, 2)
        # @test approxeq(inv(it)(south_utm), (inv(it)(south_utm)[1], inv(it)(south_utm)[2]))
    end
    
    @testset "UnitSphereFromUTM and UTMFromUnitSphere" begin
        us = SST.UnitSphereFromUTM()
        utm = SST.UTMFromUnitSphere()
        
        utm_coord = (681231.0, 5.643e6, 32, 1)
        sphere = us(utm_coord)
        @test sphere isa UnitSphericalPoint
        @test approxeq(utm(sphere), utm_coord)
        
        @test inv(us) === utm
        @test inv(utm) === us
    end
    
    @testset "with_transform" begin
        tree = SST.UTMTree(UTM_TEST_X, UTM_TEST_Y)
        
        result = SST.with_transform(tree) do t
            SST.rootnode(t)
        end
        @test result isa SST.TreeNode
        
        result2 = SST.with_transform(tree.trans) do t
            t((681231.0, 5.643e6, 32, 1))
        end
        @test result2 isa UnitSphericalPoint
    end
end
