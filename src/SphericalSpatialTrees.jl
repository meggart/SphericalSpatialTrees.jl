module SphericalSpatialTrees

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

function find_nearest(tree, point)
    cur = Ref((Inf, -1))
    point3 = UnitSphereFromGeographic()(point)
    pred = p -> begin
        _contains(p, point3)
    end
    do_query(pred, rootnode(tree)) do i
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
