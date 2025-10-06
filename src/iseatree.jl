using GeometryOps.UnitSpherical: SphericalCap, UnitSphereFromGeographic, GeographicFromUnitSphere, _intersects, _merge
using .NativeISEA: ISEA10, ISEA, RotateISEA, InvRotateISEA, PickPlane
using StaticArrays: @SVector
import YAXArrayBase
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
get_projection(t::ISEACircleTree) = inv(RotateISEA() ∘ ISEA10(t.isea))
nchild(n::ISEACircleTree) = 10
getchild(n::ISEACircleTree) = (getchild(n,i) for i in 1:nchild(n))
rootnode(t::ISEACircleTree) = t
extent(::ISEACircleTree) = SphericalCap(UnitSphereFromGeographic()((0.0,0.0)),π)
isleaf(::ISEACircleTree) = false
node_extent(t::ISEACircleTree) = extent(t)
nleaf(t::ISEACircleTree) = 10*2^(2*t.resolution)

function Base.show(io::IO, tree::ISEACircleTree)
    elements = 10 * 4^tree.resolution
    print(io, "$(elements)-leaf ISEACircleTree($(tree.resolution))")
end

function get_gridextent(tree::ISEACircleTree, xr::AbstractUnitRange, yr::AbstractUnitRange, nr::AbstractUnitRange)
    mapreduce(_merge,nr) do n
        c = getchild(tree, n)
        t = TreeNode(c.grid, TreeIndex((first(xr), last(xr) + 1), (first(yr), last(yr) + 1)))
        node_extent(t)
    end
end

struct ISEATag
    resolution::Int
    plane::Int
end
function linind(tag::ISEATag, t::TreeNode)
    LinearIndices((2^tag.resolution, 2^tag.resolution, 10))[t.index.x[1], t.index.y[1], tag.plane]
end

function get_xyranges(t::ISEACircleTree)
    t.xr,t.yr
end


getchild(t::ISEACircleTree, i::AbstractVector) = getchild(t,only(i))
function getchild(t::ISEACircleTree, i)
    xr,yr = get_xyranges(t)
    t1 = inv(ISEA10(t.isea)) ∘ InvRotateISEA() ∘ PickPlane(i)
    rootnode(RegularGridTree(xr,yr,t1,ISEATag(t.resolution,i)))
end

function index_to_cartesian(i::Integer, t::ISEACircleTree)
    n = 2^t.resolution
    CartesianIndices((n,n,10))[i].I
end
function index_to_unitsphere(i::Integer,t::ISEACircleTree)
    xr,yr = get_xyranges(t)
    halfstep = step(xr) / 2
    i,j,k = index_to_cartesian(i,t)
    x = xr[i] + halfstep
    y = yr[j] + halfstep
    getchild(t,k).grid.trans((x, y))
end

function index_to_polygon_unitsphere(i::Integer,t::ISEACircleTree)
    i,j,k = index_to_cartesian(i,t)
    index_to_polygon_unitsphere(CartesianIndex(i,j,k),t)
end

function index_to_polygon_unitsphere(i::CartesianIndex,t::ISEACircleTree)
    xr,yr = get_xyranges(t)
    i,j,k = i.I
    poly = #=@SVector =#[(xr[i], yr[j]), (xr[i+1], yr[j]), (xr[i+1], yr[j+1]), (xr[i], yr[j+1]), (xr[i], yr[j])]
    # Get the child face at index `k`, 
    # that has a transformation back to the unit sphere.
    getchild(t,k).grid.trans.(poly)
end

"""
    get_subtree(source_tree, target_chunk, target_tree)

Expands the target chunk to to the full subtree containing all grid cells at the target resolution.
"""
function get_subtree(tree::ISEACircleTree, target_indices)
    ix,iy,n = target_indices
    if length(n) > 1
        return rootnode(tree)
    end
    gridtree = getchild(tree, n).grid
    xsub = gridtree.x[ix]
    ysub = gridtree.y[iy]
    trsmall = RegularGridTree(xsub,ysub,gridtree.trans,gridtree.tag)
    rootnode(trsmall)
end


function ProjectionTarget(::Type{ISEACircleTree},target_resolution, chunk_resolution;iseaargs=())
    isea = ISEA(iseaargs...)
    tree = ISEACircleTree(isea,target_resolution)
    chunktree = ISEACircleTree(chunk_resolution)
    ProjectionTarget(tree,chunktree)
end




function create_dataset(target::ProjectionTarget{<:ISEACircleTree}, 
    path; arrayname=:layer, arraymeta=Dict(), datasetmeta=Dict(), backend=:zarr, output_datatype=Float64, kwargs...)

    back = YAXArrayBase.backendfrompath(path;driver=backend)
    axnames = ("dggs_i","dggs_j","dggs_n")
    cs = 2^(target.tree.resolution-target.chunktree.resolution)
    group = YAXArrayBase.create_dataset(
            back,
            path,
            datasetmeta,
            axnames,
            (0:(2^target.tree.resolution-1),0:(2^target.tree.resolution-1),0:9),
            (Dict(),Dict(),Dict()),
            (output_datatype,),
            (string(arrayname),),
            (axnames,),
            (arraymeta,),
            ((cs,cs,1),);
            kwargs...
        )
    group[string(arrayname)]
end

function ProjectionSource(::Type{<:ISEACircleTree}, ar, spatial_dims = (:dggs_i,:dggs_j,:dggs_n))
    isea = ISEA()
    nx,ny,n = size(ar)
    @assert n == 10 "The target must have 10 faces, got $n"
    @assert nx == ny "All ISEA faces should be square"
    flev = log2(nx)
    @assert isinteger(lev) "Size of face square should must have power of 2 length"
    lev = Int(flev)
    tree = ISEACircleTree(isea,lev)
    chunks = map(eachchunk(ar.data).chunks,DD.dims(ar)) do c,d
        DD.rebuild(d,c)
    end
    hsx = step(tree.xr)/2
    hsy = step(tree.yr)/2
    xrmid = range(first(tree.xr)+hsx,last(tree.xr)-hsx,nx)
    yrmid = range(first(tree.yr)+hsy,last(tree.yr)-hsy,ny)
    lookups = DD.Dim{:dggs_i}(xrmid), DD.Dim{:dggs_j}(yrmid),  DD.Dim{:dggs_n}(1:10)
    lookups = DD.format.(lookups)
    xchunks,ychunks,nchunks = DD.dims(chunks,spatial_dims)
    chunkres = Int(log2(length(xchunks)))
    chunktree = ISEACircleTree(isea,chunkres)
    ProjectionSource(ar,tree,chunktree,lookups,chunks)
end

function indices_from_chunk(s::ProjectionSource{<:Any,<:ISEACircleTree}, target_chunk)
    inds = index_to_cartesian(target_chunk,s.chunktree)
    chunkrange = map(getindex,s.chunks,inds)
    map(chunkrange) do cr
        Colon()(extrema(cr)...)
    end
end

# """
#     indices_from_chunk(s::ISEASourceOrTarget, target_indices)

# For a given index from the chunk tree, returns the cartesian index ranges 
# in the high-resolution tree
# """
# function indices_from_chunk(s::ISEASource, target_chunk)
#     ix, iy, n = lin_to_cart(target_chunk + 1, s.chunktree)
#     resolution_difference = s.tree.resolution - s.chunktree.resolution
#     fac = 2^resolution_difference
#     ix1 = (ix - 1) * fac + 1
#     iy1 = (iy - 1) * fac + 1
#     ix2 = ix1 + fac
#     iy2 = iy1 + fac
#     return ix1:ix2, iy1:iy2, n
# end