module SphericalSpatialTrees

import GeoFormatTypes as GFT, GeometryOps as GO, GeoInterface as GI

export RegularGridTree


#Native implementation of isea projection
include("nativeisea.jl")
#Implementation of a spherical tree in regular grids
include("RegularGridTree.jl")
#Implementation of a spherical tree for isea grid on 10 diamonds
include("iseatree.jl")
#Reprojection code
include("LazyProjection.jl")



function index_to_lonlat(i::Integer, t)
    uind = index_to_unitsphere(i, t)
    GeographicFromUnitSphere()(uind)
end

function index_to_polygon_lonlat(i, t)
    unitsphere_poly = index_to_polygon_unitsphere(i, t)
    lonlat = GeographicFromUnitSphere().(unitsphere_poly)
    lonlat_poly = GI.Polygon(#=@SVector =# [GI.LineString(lonlat; crs = GFT.EPSG(4326))]; crs = GFT.EPSG(4326))
    return lonlat_poly
end

function find_nearest(tree, point)
    cur = Ref((Inf, -1))
    point3 = UnitSphereFromGeographic()(point)
    pred = p -> begin
        _contains(p, point3)
    end
    depth_first_search(pred, rootnode(tree)) do i
        cur_min = first(cur[])
        dist = norm(point3 - index_to_unitsphere(i, tree))
        if dist < cur_min
            cur[] = (dist, i)
        end
    end
    cur[]
end


sanitize_predicate(pred::SphericalCap) = Base.Fix1(_intersects, pred)

end