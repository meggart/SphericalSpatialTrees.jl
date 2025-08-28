using SphericalSpatialTrees: RegularGridTree, TreeNode, nchild, getchild, rootnode, extent, isleaf, node_extent, 
    linind
using GeometryOps.UnitSpherical: SphericalCap, UnitSphereFromGeographic, GeographicFromUnitSphere, _intersects
using .NativeISEA: ISEA10, ISEA, InvISEA10
using StaticArrays: @SVector

struct ISEACircleTree
    isea::ISEA
    resolution::Int
end
ISEACircleTree(resolution::Int) = ISEACircleTree(ISEA(), resolution)
nchild(n::ISEACircleTree) = 10
getchild(n::ISEACircleTree) = (getchild(n,i) for i in 1:nchild(n))
rootnode(t::ISEACircleTree) = t
extent(::ISEACircleTree) = SphericalCap(UnitSphereFromGeographic()((0.0,0.0)),π)
isleaf(::ISEACircleTree) = false
node_extent(t::ISEACircleTree) = extent(t)
struct _RotatedISEAtoUnitSphere{T}
    i::Int
    t::T
end
function (r::_RotatedISEAtoUnitSphere)((x,y))
    x, y = NativeISEA.itrans_matrix * @SVector([x, y * sqrt(3) / 2])
    r.t((r.i, x, y))
end
_RotatedISEAtoUnitSphere(i::Int) = _RotatedISEAtoUnitSphere(i, inv(ISEA10()))



struct _RotatedISEAtolatlon{T}
    i::Int
    t::T
end
function (r::_RotatedISEAtolatlon)((x,y))
    x, y = NativeISEA.itrans_matrix * @SVector([x, y * sqrt(3) / 2])
    res = r.t((r.i, x, y))
    GeographicFromUnitSphere()(res)
end
_RotatedISEAtolatlon(i::Int) = _RotatedISEAtolatlon(i, inv(ISEA10()))
struct ISEATag
    resolution::Int
end
function linind(tag::ISEATag, t::TreeNode)
    #@show t.index.x[1]-1,t.index.y[1]-1,t.grid.trans.i-1,tag.resolution
    # @assert t.index.x[1]==t.index.x[2]-1
    # @assert t.index.y[1]==t.index.y[2]-1
    LinearIndices((2^tag.resolution,2^tag.resolution,10))[t.index.x[1], t.index.y[1], t.grid.trans.i] - 1
    #Int(Cell10(t.index.x[1] - 1, t.index.y[1] - 1, t.grid.trans.i - 1, tag.resolution))
    #t.index.x[1] + (t.index.y[1]-1)*(length(t.grid.x)-1)
end

function getchild(t::ISEACircleTree, i)
    n = 2^t.resolution
    xr = range(0, 1, length=n + 1)
    yr = range(0,1,length=n+1)
    t1 = _RotatedISEAtoUnitSphere(i,InvISEA10(t.isea))
    rootnode(RegularGridTree(xr,yr,t1,ISEATag(t.resolution)))
end

function index_to_unitsphere(i::Integer,t::ISEACircleTree)
    n = 2^t.resolution
    xr = range(0, 1, length=n + 1)
    yr = range(0, 1, length=n + 1)
    halfstep = step(xr) / 2
    i,j,n = CartesianIndices((n,n,10))[i+1].I
    x = xr[i+1] + halfstep
    y = yr[j+1] + halfstep
    getchild(t,n).grid.trans((x, y))
end


function index_to_lonlat(i::Integer,t)
    uind = index_to_unitsphere(i,t)
    GeographicFromUnitSphere()(uind)
end
