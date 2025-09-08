using GeometryOps.UnitSpherical: SphericalCap, UnitSphereFromGeographic, GeographicFromUnitSphere, _intersects
using .NativeISEA: ISEA10, ISEA, RotateISEA, InvRotateISEA, PickPlane
using StaticArrays: @SVector
using YAXArrays: YAXArray, setchunks, Dataset, savedataset
using FillArrays: Fill

struct ISEACircleTree
    isea::ISEA{Float64}
    resolution::Int
    xr::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    yr::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
end
ISEACircleTree(isea::ISEA,resolution::Int) = ISEACircleTree(isea,resolution, range(0.0,1.0,length=2^resolution+1),range(0.0,1.0,length=2^resolution+1))
ISEACircleTree(resolution::Int) = ISEACircleTree(ISEA(),resolution)
gridsize(t::ISEACircleTree) = (2^t.resolution, 2^t.resolution, 10)
Base.ndims(t::ISEACircleTree) = 3
get_projection(t::ISEACircleTree) = RotateISEA() ∘ ISEA10(t.isea)
nchild(n::ISEACircleTree) = 10
getchild(n::ISEACircleTree) = (getchild(n,i) for i in 1:nchild(n))
rootnode(t::ISEACircleTree) = t
extent(::ISEACircleTree) = SphericalCap(UnitSphereFromGeographic()((0.0,0.0)),π)
isleaf(::ISEACircleTree) = false
node_extent(t::ISEACircleTree) = extent(t)
nleaf(t::ISEACircleTree) = 10*2^(2*t.resolution)


# struct _RotatedISEAtoUnitSphere{T}
#     i::Int
#     t::T
# end
# function (r::_RotatedISEAtoUnitSphere)((x,y))
#     #x, y = NativeISEA.itrans_matrix * @SVector([x, y * sqrt(3) / 2])
#     r.t((r.i, x, y))
# end
# # _RotatedISEAtoUnitSphere(i::Int) = _RotatedISEAtoUnitSphere(i, inv(ISEA10()))


# struct _RotatedISEAtolatlon{T}
#     i::Int
#     t::T
# end
# function (r::_RotatedISEAtolatlon)((x,y))
#     #x, y = NativeISEA.itrans_matrix * @SVector([x, y * sqrt(3) / 2])
#     res = r.t((r.i, x, y))
#     GeographicFromUnitSphere()(res)
# end
# _RotatedISEAtolatlon(i::Int) = _RotatedISEAtolatlon(i, inv(ISEA10()))
struct ISEATag
    resolution::Int
    plane::Int
end
function linind(tag::ISEATag, t::TreeNode)
    #@show t.index.x[1]-1,t.index.y[1]-1,t.grid.trans.i-1,tag.resolution
    # @assert t.index.x[1]==t.index.x[2]-1
    # @assert t.index.y[1]==t.index.y[2]-1
    LinearIndices((2^tag.resolution,2^tag.resolution,10))[t.index.x[1], t.index.y[1], tag.plane]
    #Int(Cell10(t.index.x[1] - 1, t.index.y[1] - 1, t.grid.trans.i - 1, tag.resolution))
    #t.index.x[1] + (t.index.y[1]-1)*(length(t.grid.x)-1)
end

function get_xyranges(t::ISEACircleTree)
    t.xr,t.yr
end

function getchild(t::ISEACircleTree, i)
    xr,yr = get_xyranges(t)
    t1 = inv(ISEA10(t.isea)) ∘ InvRotateISEA() ∘ PickPlane(i)
    rootnode(RegularGridTree(xr,yr,t1,ISEATag(t.resolution,i)))
end

function lin_to_cart(i::Integer,t::ISEACircleTree)
    n = 2^t.resolution
    CartesianIndices((n,n,10))[i].I
end
function index_to_unitsphere(i::Integer,t::ISEACircleTree)
    xr,yr = get_xyranges(t)
    halfstep = step(xr) / 2
    i,j,k = lin_to_cart(i,t)
    x = xr[i] + halfstep
    y = yr[j] + halfstep
    getchild(t,k).grid.trans((x, y))
end


function index_to_lonlat(i::Integer,t)
    uind = index_to_unitsphere(i,t)
    GeographicFromUnitSphere()(uind)
end

"""
    get_subtree(source_tree, target_chunk, target_tree)

Expands the target chunk to to the full subtree containing all grid cells at the target resolution.
"""
function get_subtree(source_tree::ISEACircleTree, target_chunk, target_tree::ISEACircleTree)
    ix, iy, n = lin_to_cart(target_chunk + 1, source_tree)
    resolution_difference = target_tree.resolution - source_tree.resolution
    fac = 2^resolution_difference
    ix1 = (ix - 1) * fac + 1
    iy1 = (iy - 1) * fac + 1
    ix2 = ix1 + fac
    iy2 = iy1 + fac
    gridtree = getchild(target_tree, n).grid
    xsub = gridtree.x[ix1:ix2]
    ysub = gridtree.y[iy1:iy2]
    trsmall = RegularGridTree(xsub,ysub,gridtree.trans,gridtree.tag)
    rootnode(trsmall)
end

function get_subindices(source_tree::ISEACircleTree, target_chunk, target_tree::ISEACircleTree)
    ix, iy, n = lin_to_cart(target_chunk, source_tree)
    resolution_difference = target_tree.resolution - source_tree.resolution
    fac = 2^resolution_difference
    ix1 = (ix - 1) * fac + 1
    iy1 = (iy - 1) * fac + 1
    ix2 = ix1 + fac -1 
    iy2 = iy1 + fac -1
    ix1:ix2, iy1:iy2, n
end

struct ProjectionTarget{T,CT}
    tree::T
    chunktree::CT
end
function ProjectionTarget(::Type{ISEACircleTree},target_resolution, chunk_resolution;iseaargs=())
    isea = ISEA(iseaargs...)
    tree = ISEACircleTree(isea,target_resolution)
    chunktree = ISEACircleTree(chunk_resolution)
    ProjectionTarget(tree,chunktree)
end

function create_dataset(target::ProjectionTarget{<:ISEACircleTree}, 
    path; arrayname=:layer, arraymeta=Dict(), datasetmeta=Dict())
    outdims = (DD.Dim{:dggs_i}(0:(2^target.tree.resolution-1)),
        DD.Dim{:dggs_j}(0:(2^target.tree.resolution-1)),
        DD.Dim{:dggs_n}(0:9),
    )
    a = YAXArray(outdims, Fill(0.0,length.(outdims)),arraymeta)
    cs = 2^(target.tree.resolution-target.chunktree.resolution)
    p = Symbol(arrayname)=>setchunks(a,(cs,cs,1))
    p = (p,)
    ds = Dataset(;properties=datasetmeta,p...)
    ds = savedataset(ds,path = path, skeleton=true, overwrite=true)
end