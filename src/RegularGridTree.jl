using GeometryOps.UnitSpherical: UnitSpherical, _contains, _intersects, UnitSphereFromGeographic,
    spherical_distance, SphericalCap, UnitSphericalPoint
import GeometryOps.SpatialTreeInterface: nchild, getchild, isleaf, child_indices_extents,
    query, sanitize_predicate, node_extent, depth_first_search
import GeometryOps: extent
import StaticArrays: @SVector
import GeometryOps.Extents: Extent, bounds
import LinearAlgebra: norm, cross, dot
import DimensionalData as DD


struct RegularGridTree{DX,DY,T,S}
    x::DX
    y::DY
    trans::T
    tag::S
end
Base.ndims(t::RegularGridTree) = 2
gridsize(t::RegularGridTree) = (length(t.x)-1, length(t.y)-1)
function get_gridextent(t::RegularGridTree, xr::AbstractUnitRange, yr::AbstractUnitRange)
    t = TreeNode(t, TreeIndex((first(xr), last(xr) + 1), (first(yr), last(yr) + 1)))
    node_extent(t)
end
get_projection(t::RegularGridTree) = t.trans
function DD.dims(r::RegularGridTree)
    xmid, ymid = map((r.x, r.y)) do d
        (d[1:(end-1)] .+ d[2:end]) ./ 2
    end
    DD.X(xmid), DD.Y(ymid)
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
    xr, yr = map(x->boundrangefromcenters(x.val), ar_spatial_dims)
    return RegularGridTree(xr, yr, transform)
end

get_tag(r::RegularGridTree) = r.tag
function extent(r::RegularGridTree, x1, x2, y1, y2)
    Extent(X=(r.x[x1], r.x[x2]), Y=(r.y[y1], r.y[y2]))
end
nlevel(r::RegularGridTree) = max(ceil(Int, log(length(r.x))), ceil(Int, log(length(r.y))))


abstract type AbstractTreeIndex end
struct TreeIndex <: AbstractTreeIndex
    x::Tuple{Int,Int}
    y::Tuple{Int,Int}
end
is_valid_index(t::TreeIndex) = (first(t.x) < last(t.x) && first(t.y) < last(t.y))

_isone(t) = last(t) - first(t) == 1
_nchild(i::TreeIndex) = 4 ÷ (_isone(i.x) + _isone(i.y) + 1)
function split_4(index::TreeIndex)
    l, r = split_half(index.x)
    u, d = split_half(index.y)
    TreeIndex(l, u), TreeIndex(l, d), TreeIndex(r, u), TreeIndex(r, d)
end

function split_half((x1, x2))
    middle = x1 + div(x2 - x1, 2)
    (x1, middle), (middle, x2)
end
struct TreeNode{T,I}
    grid::T
    index::I
    function TreeNode(grid, index::AbstractTreeIndex)
        if !is_valid_index(index)
            println(index)
            error()
        end
        new{typeof(grid),typeof(index)}(grid, index)
    end
end
linind(t::TreeNode) = linind(get_tag(t.grid), t)
linind(::Nothing, t::TreeNode) = linind(t.grid, t.index)
function linind(grid::RegularGridTree, index::TreeIndex)
    LinearIndices((length(grid.x)-1, length(grid.y)-1))[index.x[1], index.y[1]]
end
function correct_index(i, index)
    i = ifelse(_isone(index.x), i + 2, i)
    i = ifelse(_isone(index.y), i * 2, i)
    i
end
function getchild(t::TreeNode, i)
    i2 = correct_index(i, t.index)
    index = split_4(t.index)[i2]
    try
        TreeNode(t.grid, index)
    catch e
        @show i
        @show t.index
        @show i2
        @show split_4(t.index)
        rethrow(e)
    end
end
nchild(n::TreeNode) = _nchild(n.index)
getchild(n::TreeNode) = (getchild(n, i) for i in 1:nchild(n))
"""
    rootnode(t)

Return the root [`TreeNode`](@ref) of the spatial tree `t`, the starting point
for searches and traversals.
"""
rootnode(t::RegularGridTree) = TreeNode(t, TreeIndex((1, length(t.x)), (1, length(t.y))))
extent(node::TreeNode) = extent(node.grid, node.index)
extent(t::RegularGridTree, index::TreeIndex) = extent(t, index.x[1], index.x[2], index.y[1], index.y[2])
isleaf(node::TreeNode) = isleaf(node.index)
isleaf(index::TreeIndex) = _isone(index.x) && _isone(index.y)

node_extent(node::TreeNode) = circle_from_extent_1(extent(node), node.grid)
function child_indices_extents(tree::TreeNode)
    li = linind(tree)
    circ = circle_from_extent_1(extent(tree), tree.grid)
    ((li, circ),)
end
nleaf(t::RegularGridTree) = (length(t.x)-1)*(length(t.y)-1)

node_to_polygon_unitsphere(i::TreeNode) = node_to_polygon_unitsphere(i.grid, i.index)
function node_to_polygon_unitsphere(grid::RegularGridTree, index::TreeIndex)
    x1, x2 = index.x
    y1, y2 = index.y
    xr = grid.x
    yr = grid.y
    poly = #= @SVector =#[(xr[x1], yr[y1]), (xr[x2], yr[y1]), (xr[x2], yr[y2]), (xr[x1], yr[y2]), (xr[x1], yr[y1])]
    # Get the child face at index `k`, 
    # that has a transformation back to the unit sphere.
    grid.trans.(poly)
end

function _circle_from_pair(a, b)
    c_ = (a+b)
    c = c_/norm(c_)
    r = spherical_distance(c, a)
    SphericalCap(c, r*1.0001)
end
function _circle_from_3(a, b, c)
    n_ = cross((b-a), (c-a))
    sum(n_)==0.0 && return SphericalCap(a, 0.0)
    n = UnitSphericalPoint(n_/norm(n_))
    if dot(n,a) < 0.0
        n = -n
    end
    d = spherical_distance(a, n)
    SphericalCap(n, d*1.0001)
end

circle_from_extent_1(ex, grid) = _circle_from_extent(ex, grid.trans)

function _circle_from_extent(ex, trans)
    (x1, x2), (y1, y2) = bounds(ex)
    cx, cy = (x2 + x1) / 2, (y2 + y1) / 2
    a,b,c,d,e,f,g,h = map(trans, ((x1, y1), (x2, y1), (x2, y2), (x1, y2),(cx, y1), (x2, cy), (cx, y2), (x1, cy)))
    # Determine if we might run into 
    has_large = any(((a,e),(b,f),(c,g),(d,h))) do (x,y)
        ang = dot(x,y)
        ang < 0.7 || ang==1.0
    end
    cap = if has_large
        z = trans((cx,cy))
        alld = map(p->spherical_distance(z, p), (a,b,c,d,e,f,g,h))
        r = reduce(max, alld)
        SphericalCap(z, r)
    else
        for points in ((a, c, b, d), (b, d, a, c), (a, b, c, d), (b, c, a, d), (c, d, a, b), (d, a, b, c))
            cap = _circle_from_pair(points[1], points[2])
            if _contains(cap, points[3]) && _contains(cap, points[4])
                return cap
            end
        end
        for points in ((a, b, c, d), (b, c, d, a), (c, d, a, b), (d, a, b, c))
            cap = _circle_from_3(points[1], points[2], points[3])
            if _contains(cap, points[4])
                return cap
            end
        end
        mapreduce(_merge,(a,c,b,d)) do p
            SphericalCap(p,0.0)
        end
    end
    #The following is done to not miss intersections through numerical inaccuracies
    SphericalCap(cap.point, cap.radius*1.0001)
end

get_step(x::AbstractRange) = step(x)
get_step(x) = length(x) > 1 ? (x[2] - x[1]) : one(eltype(x))

"""
    index_to_native_coords(i, t)

Return the coordinates of the center of cell index `i` in the native (unprojected)
coordinate system of tree `t`.
"""
function index_to_native_coords(i, t::RegularGridTree)
    xhalfstep = get_step(t.x) / 2
    yhalfstep = get_step(t.y) / 2
    i, j = index_to_cartesian(i, t)
    x = t.x[i] + xhalfstep
    y = t.y[j] + yhalfstep
    x, y
end

function index_to_lonlat(i::Integer, t::RegularGridTree{<:Any,<:Any,UnitSphereFromGeographic})
    index_to_native_coords(i, t)
end

"""
    index_to_cartesian(i, t)

Convert a linear index `i` into the Cartesian index (tuple of subscripts) of the
corresponding cell of tree `t`.
"""
index_to_cartesian(i::Integer, t::RegularGridTree) = CartesianIndices((length(t.x) - 1, length(t.y) - 1))[i].I

"""
    index_to_unitsphere(i, t, projfunc=get_projection(t))

Return the position of the center of cell index `i` of tree `t` as a
`UnitSphericalPoint` on the unit sphere.
"""
function index_to_unitsphere(i::Integer, t, projfunc=get_projection(t))
    coords = index_to_native_coords(i, t)
    projfunc(coords)
end



function TreeNode(tree::RegularGridTree, targetinds::Tuple)
    r1, r2 = targetinds
    ix1, ix2 = first(r1), last(r1)
    iy1, iy2 = first(r2), last(r2)
    TreeNode(tree, TreeIndex((ix1, ix2+1), (iy1, iy2+1)))
end

"""
    ProjectionSource(::Type{<:RegularGridTree}, ar, spatial_dims=(DD.XDim, DD.YDim))

Create a regridding source from a `DD.AbstractDimArray` with regular grid
dimensions (by default `X` and `Y`). The source chunk tree is derived from the
chunking of `ar`.
"""
function ProjectionSource(::Type{<:RegularGridTree}, ar, spatial_dims=(DD.XDim, DD.YDim))
    tree = RegularGridTree(ar, spatial_dims)
    lookups = map(DD.format, DD.dims(ar, spatial_dims))
    chunks = map(eachchunk(ar.data).chunks, DD.dims(ar)) do c, d
        DD.rebuild(d, c)
    end
    xchunks, ychunks = DD.dims(chunks, spatial_dims)
    xchunkbnds = vcat(tree.x[first.(xchunks.val)], last(tree.x))
    ychunkbnds = vcat(tree.y[first.(ychunks.val)], last(tree.y))
    chunktree = RegularGridTree(xchunkbnds, ychunkbnds)
    ProjectionSource(ar, tree, chunktree, lookups, chunks)
end

#Compute indices given a chunk index for the high-resolution tree
function indices_from_chunk(s::ProjectionSource{<:Any,<:RegularGridTree}, target_chunk)
    inds = index_to_cartesian(target_chunk, s.chunktree)
    chunkrange = map(getindex, s.chunks, inds)
    map(chunkrange) do cr
        Colon()(extrema(cr)...)
    end
end


"""
    ProjectionTarget(::Type{<:RegularGridTree}, x, y, trans=UnitSphereFromGeographic(); chunksize=256)

Create a regridding target on a regular grid with cell boundaries `x` and `y`.
`chunksize` sets the number of cells per chunk in each direction.
"""
function ProjectionTarget(::Type{<:RegularGridTree}, x, y, trans=UnitSphereFromGeographic(); chunksize=256)
    tree = RegularGridTree(x, y, trans)
    xchunk = x[1:chunksize:end]
    ychunk = y[1:chunksize:end]
    chunktree = RegularGridTree(xchunk, ychunk, trans)
    ProjectionTarget(tree, chunktree)
end
Base.ndims(::Type{RegularGridTree}) = 2
