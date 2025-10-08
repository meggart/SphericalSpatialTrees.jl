using GeometryOps.UnitSpherical: UnitSpherical, _contains, _intersects, UnitSphereFromGeographic,
    spherical_distance, SphericalCap, UnitSphericalPoint
import GeometryOps.SpatialTreeInterface: nchild, getchild, isleaf, child_indices_extents,
    query, sanitize_predicate, node_extent, depth_first_search
import GeometryOps: extent
import StaticArrays: @SVector
import GeometryOps.Extents: Extent, bounds
import LinearAlgebra: norm
import DimensionalData as DD


struct RegularGridTree{DX,DY,T,S}
    x::DX
    y::DY
    trans::T
    tag::S
end
Base.ndims(t::RegularGridTree) = 2
gridsize(t::RegularGridTree) = (length(t.x)-1,length(t.y)-1)
function get_gridextent(t::RegularGridTree, xr::AbstractUnitRange, yr::AbstractUnitRange)
    t = TreeNode(t, TreeIndex((first(xr), last(xr) + 1), (first(yr), last(yr) + 1)))
    node_extent(t)
end
get_projection(t::RegularGridTree) = t.trans
function DD.dims(r::RegularGridTree)
    xmid,ymid = map((r.x,r.y)) do d
        (d[1:end-1] .+ d[2:end])./2
    end
    DD.X(xmid),DD.Y(ymid)
end

function Base.show(io::IO, tree::RegularGridTree)
    # Check if this is a compact display (when used within other show methods)
    compact = get(io, :compact, false)
    
    if compact
        # Compact format for use within other show methods
        print(io, "RegularGridTree($(length(tree.x))×$(length(tree.y)))")
    else
        # Full format with copy-pastable constructor on first line
        print(io, "RegularGridTree($(length(tree.x))×$(length(tree.y)) array, chunksize)")
        
        # Additional lines in cyan color with supplementary information
        if get(io, :color, false)
            print(io, "\n\e[36m")  # Start cyan color
            print(io, "dimensions: $(length(tree.x)-1)×$(length(tree.y)-1)")
            print(io, "\e[0m")     # Reset color
        else
            print(io, "\ndimensions: $(length(tree.x)-1)×$(length(tree.y)-1)")
        end
    end
end
"""
    RegularGridTree(x, y, transform=UnitSphereFromGeographic())

Constructs a simple 2D rectangular grid as a quad tree adhering to the SpatialTreeInterface. `x` and `y` represent 
the boundaries of the grid cells, i.e. `length(x) == n_grid_x + 1`. 
Any transformation can be passed that converts static vector values of [x,y] to a `UnitSphericalPoint`.
"""
RegularGridTree(x, y, transform=UnitSphereFromGeographic()) = RegularGridTree(x, y, transform, nothing)

"""
    RegularGridTree(ar::DD.AbstractDimArray,spatial_dims;transform=UnitSphereFromGeographic())

Convenience constructor to create a unitspherical spatial search tree from an AbstractDimArray. The spatial dimension
names can be passed as a tuple of symbols. Tries to guess cell boundaries. 
"""
function RegularGridTree(ar::DD.AbstractDimArray, spatial_dims=(DD.XDim, DD.YDim); transform=UnitSphereFromGeographic())
    ar_spatial_dims = DD.dims(ar, spatial_dims)
    if isnothing(ar_spatial_dims) || any(isnothing, ar_spatial_dims)
        dimstrings = map(spatial_dims) do d
            dtype = if d isa Type
                d
            else
                typeof(d)
            end
            sprint(io -> Base.show_type_name(io, Core.typename(dtype)))
        end

        error("""
            You requested the spatial dims
            `$(join(dimstrings, ", "))`
            but they could not be found in your dimarray with dims
            `$(join(map(DD.name, DD.dims(ar)), ", "))`.
            Please pass `spatial_dims` that exist within the array.
        """)
    end
    xr, yr = map(boundrangefromcenters, ar_spatial_dims)
    return RegularGridTree(xr, yr, transform)
end

get_tag(r::RegularGridTree) = r.tag
function extent(r::RegularGridTree, x1, x2, y1, y2)
    Extent(X=(r.x[x1], r.x[x2]), Y=(r.y[y1], r.y[y2]))
end
nlevel(r::RegularGridTree) = max(ceil(Int, log(length(r.x))), ceil(Int, log(length(r.y))))

struct TreeIndex
    x::Tuple{Int,Int}
    y::Tuple{Int,Int}
    function TreeIndex(x::Tuple{Int,Int}, y::Tuple{Int,Int})
        (last(x) <= first(x) || last(y) <= first(y)) && error()
        new(x, y)
    end
end

_isone(t) = last(t) - first(t) == 1
_nchild(i::TreeIndex) = 4 ÷ (_isone(i.x) + _isone(i.y) + 1)
split_left((x1, x2)) = (x1, x1 + div(x2 - x1, 2))
split_right((x1, x2)) = (x1 + div(x2 - x1, 2), x2)
struct TreeNode{T<:RegularGridTree}
    grid::T
    index::TreeIndex
end
linind(t::TreeNode) = linind(get_tag(t.grid), t)
function linind(::Nothing, t::TreeNode)
    LinearIndices((length(t.grid.x)-1, length(t.grid.y)-1))[t.index.x[1], t.index.y[1]]
end
function correct_index(i, index)
    i = ifelse(_isone(index.x), i + 2, i)
    i = ifelse(_isone(index.y), i * 2, i)
    i
end
function getchild(t::TreeNode, i)
    i = correct_index(i, t.index)
    i == 1 ? TreeNode(t.grid, TreeIndex(split_left(t.index.x), split_left(t.index.y))) :
    i == 2 ? TreeNode(t.grid, TreeIndex(split_left(t.index.x), split_right(t.index.y))) :
    i == 3 ? TreeNode(t.grid, TreeIndex(split_right(t.index.x), split_left(t.index.y))) :
    i == 4 ? TreeNode(t.grid, TreeIndex(split_right(t.index.x), split_right(t.index.y))) :
    throw(ArgumentError("TreeNodes can only be indexed with indices 1:4"))
end
nchild(n::TreeNode) = _nchild(n.index)
getchild(n::TreeNode) = (getchild(n, i) for i in 1:nchild(n))
rootnode(t::RegularGridTree) = TreeNode(t, TreeIndex((1, length(t.x)), (1, length(t.y))))
extent(node::TreeNode) = extent(node.grid, node.index.x[1], node.index.x[2], node.index.y[1], node.index.y[2])
isleaf(node::TreeNode) = _isone(node.index.x) && _isone(node.index.y)
node_extent(node::TreeNode) = circle_from_extent_1(extent(node), node.grid.trans)
function child_indices_extents(tree::TreeNode) 
    li = linind(tree)
    circ = circle_from_extent_1(extent(tree), tree.grid.trans)
    ((li, circ) ,)
end
nleaf(t::RegularGridTree) = (length(t.x)-1)*(length(t.y)-1)

function node_to_polygon_unitsphere(i::TreeNode)
    x1,x2 = i.index.x
    y1,y2 = i.index.y
    xr = i.grid.x
    yr = i.grid.y
    poly = #= @SVector =# [(xr[x1], yr[y1]), (xr[x2], yr[y1]), (xr[x2], yr[y2]), (xr[x1], yr[y2]), (xr[x1], yr[y1])]
    # Get the child face at index `k`, 
    # that has a transformation back to the unit sphere.
    i.grid.trans.(poly)
end

function circle_from_extent_1(ex, trans)
    (x1, x2), (y1, y2) = bounds(ex)
    cx,cy = (x2 + x1) / 2, (y2 + y1) / 2
    a, b, c, d, e, f, g, h = map(trans, ((x1, y1), (x2, y1), (x2, y2), (x1, y2),(cx,y1),(x2,cy),(cx,y2),(x1,cy)))
    z = trans((cx,cy))
    alld = map(p->spherical_distance(z,p), (a, b, c, d,e,f,g,h))
    r = reduce(max, alld)
    #The following is done to not miss intersections through numerical inaccuracies
    res = SphericalCap(z, r*1.0001)
    # if !all(_contains.((res,), (a,b,c,d)))
    #     @show a,b,c,d,e,f,g,h
    #     @show e
    #     @show alld
    #     error()
    # end
    res
end

get_step(x::AbstractRange) = step(x)
get_step(x) = length(x) > 1 ? (x[2] - x[1]) : one(eltype(x))

function index_to_lonlat(i::Integer,t::RegularGridTree{<:Any,<:Any,UnitSphereFromGeographic})
    xhalfstep = get_step(t.x) / 2
    yhalfstep = get_step(t.y) / 2
    i, j = index_to_cartesian(i, t)
    x = t.x[i] + xhalfstep
    y = t.y[j] + yhalfstep
    x,y
end

index_to_cartesian(i::Integer, t::RegularGridTree) = CartesianIndices((length(t.x) - 1, length(t.y) - 1))[i].I

function index_to_unitsphere(i::Integer,t::RegularGridTree)
    xhalfstep = step(t.x) / 2
    yhalfstep = step(t.y) / 2
    i, j = index_to_cartesian(i, t)
    x = t.x[i] + xhalfstep
    y = t.y[j] + yhalfstep
    t.trans((x,y))
end
# function circle_from_extent_2(ex,trans)
#     (x1, x2), (y1, y2) = bounds(ex)
#     points = map(trans, ((x1, y1), (x2, y1), (x2, y2), (x1, y2)))
#     cap1 = SphericalCap(points[1], points[2], points[3])
#     cap2 = SphericalCap(points[2], points[3], points[4])
#     UnitSpherical._merge(cap1, cap2)
# end


function TreeNode(tree::RegularGridTree,targetinds)
    r1,r2 = targetinds
    ix1,ix2 = first(r1), last(r1)
    iy1,iy2 = first(r2), last(r2)
    TreeNode(tree, TreeIndex((ix1, ix2+1), (iy1, iy2+1)))
end

function ProjectionSource(::Type{<:RegularGridTree}, ar, spatial_dims = (DD.XDim,DD.YDim))
    tree = RegularGridTree(ar,spatial_dims)
    lookups = map(DD.format,DD.dims(ar,spatial_dims))
    chunks = map(eachchunk(ar.data).chunks,DD.dims(ar)) do c,d
        DD.rebuild(d,c)
    end
    xchunks,ychunks = DD.dims(chunks,spatial_dims)
    xchunkbnds = vcat(tree.x[first.(xchunks.val)], last(tree.x))
    ychunkbnds = vcat(tree.y[first.(ychunks.val)], last(tree.y))
    chunktree = RegularGridTree(xchunkbnds,ychunkbnds)
    ProjectionSource(ar,tree,chunktree,lookups,chunks)
end

#Compute indices given a chunk index for the high-resolution tree
function indices_from_chunk(s::ProjectionSource{<:Any,<:RegularGridTree}, target_chunk)
    inds = index_to_cartesian(target_chunk, s.chunktree)
    chunkrange = map(getindex,s.chunks,inds)
    map(chunkrange) do cr
        Colon()(extrema(cr)...)
    end
end


function ProjectionTarget(::Type{<:RegularGridTree},x,y,trans=UnitSphereFromGeographic();chunksize=256)
    tree = RegularGridTree(x,y,trans)
    xchunk = x[1:chunksize:end]
    ychunk = y[1:chunksize:end]
    chunktree = RegularGridTree(xchunk,ychunk,trans)
    ProjectionTarget(tree,chunktree)
end
Base.ndims(::Type{RegularGridTree}) = 2
