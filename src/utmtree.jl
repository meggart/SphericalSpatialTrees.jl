abstract type LazyProjTransform <: Transformation end
function (t::LazyProjTransform)(x)
    with_transform(t) do trans
        trans(x)
    end
end
allow_threaded_transformation(::LazyProjTransform) = false

create_transform(t::LazyProjTransform, ctx) = if isinv(t) 
    GeographicFromUnitSphere() ∘ Proj.Transformation("OGC:84",code(t),ctx=ctx) 
else 
    Proj.Transformation(code(t),"OGC:84",ctx=ctx) ∘ UnitSphereFromGeographic()
end
function with_transform(f,t::LazyProjTransform)
    ctx = Proj.proj_context_create()
    try 
        _tt = create_transform(t,ctx)
        f(_tt)
    finally
        Proj.proj_context_destroy(ctx)
    end
end

struct UTMTransform{T} <: Transformation
    projs::T
    ctx::Ptr{Nothing}
end
function (t::UTMTransform)((lon,lat))
    hemi = 1 + (lat < 0.0)
    zone = Int((mod(lon+180,360)) ÷ 6) + 1
    if ismissing(t.projs[zone,hemi])
        pstr = hemi==2 ? "+proj=utm +zone=$zone +south" : "+proj=utm +zone=$zone"
        t.projs[zone,hemi] = Proj.Transformation("OGC:84",pstr,ctx=t.ctx) 
    end
    x,y = t.projs[zone,hemi]((lon,lat))
    return (x,y,zone,hemi)
end
struct IUTMTransform{T} <: Transformation
    projs::T
    ctx::Ptr{Nothing}
end
function (t::IUTMTransform)((x,y,zone,hemi))
    if ismissing(t.projs[zone,hemi])
        pstr = hemi==2 ? "+proj=utm +zone=$zone +south" : "+proj=utm +zone=$zone"
        t.projs[zone,hemi] = Proj.Transformation(pstr,"OGC:84";ctx=t.ctx) 
    end
    lon,lat = t.projs[zone,hemi]((x,y))
    return (lon,lat)
end
Base.inv(t::UTMTransform) = IUTMTransform(ctx = t.ctx)
Base.inv(t::IUTMTransform) = UTMTransform(ctx = t.ctx)


function IUTMTransform(;ctx=C_NULL)
    projs = Union{Missing,Proj.Transformation}[missing for i in 1:60, j in 1:2]
    IUTMTransform(projs,ctx)
end
function UTMTransform(;ctx=C_NULL)
    projs = Union{Missing,Proj.Transformation}[missing for i in 1:60, j in 1:2]
    UTMTransform(projs,ctx)
end

struct UnitSphereFromUTM <: LazyProjTransform end
struct UTMFromUnitSphere <: LazyProjTransform end
create_transform(::UnitSphereFromUTM,ctx) = UnitSphereFromGeographic() ∘ IUTMTransform(;ctx)
create_transform(::UTMFromUnitSphere,ctx) = UTMTransform(;ctx) ∘ GeographicFromUnitSphere()
Base.inv(::UnitSphereFromUTM) = UTMFromUnitSphere()
Base.inv(::UTMFromUnitSphere) = UnitSphereFromUTM()


struct UTMIndex <: AbstractTreeIndex
    x::Tuple{Int,Int}
    y::Tuple{Int,Int}
    zone::Tuple{Int,Int}
    hemi::UInt8 # 1: North, 2: South, 3: Both
end

is_valid_index(i::UTMIndex) = last(i.x)>first(i.x) && last(i.y)>first(i.y) && last(i.zone)>first(i.zone) && (1<=i.hemi<=3)
_nchild(i::UTMIndex) = 4 ÷ (_isone(i.x) + _isone(i.y) + 1)
function split_utm_updown(hemi, y)
    if hemi == 3
        (1, y),(2,y)
    else
        y1,y2 = split_half(y)
        (hemi,y1),(hemi,y2)
    end
end
function split_utm_leftright(zone,x)
    if _isone(zone)
        x1,x2 = split_half(x)
        (zone,x1),(zone,x2)
    else
        zone1,zone2 = split_half(zone)
        (zone1,x),(zone2,x)
    end
end
function split_4(index::UTMIndex)
    ((hemi1,y1),(hemi2,y2)) = split_utm_updown(index.hemi,index.y)
    ((zone1,x1),(zone2,x2)) = split_utm_leftright(index.zone,index.x)
    UTMIndex(x1,y1,zone1,hemi1), UTMIndex(x1,y2,zone1,hemi2), UTMIndex(x2,y1,zone2,hemi1), UTMIndex(x2,y2,zone2,hemi2)    
end

struct UTMTree{DX,DY,T}
    x::DX
    y::DY
    trans::T
end
function with_transform(f::F, tree::UTMTree) where F
    with_transform(tree.trans) do trans2
        tree2 = UTMTree(tree.x,tree.y,trans2)
        f(tree2)
    end
end

Base.ndims(t::UTMTree) = 4
gridsize(t::UTMTree) = (length(t.x)-1, length(t.y)-1, 60, 2)
function get_gridextent(t::UTMTree, xr::AbstractUnitRange, yr::AbstractUnitRange, zone, hemi)
    t = TreeNode(t, UTMIndex((first(xr), last(xr) + 1), (first(yr), last(yr) + 1), (first(zone), last(zone)+1), length(hemi)==2 ? 3 : first(hemi)))
    node_extent(t)
end
get_projection(t::UTMTree) = t.trans
function DD.dims(r::UTMTree)
    xmid, ymid = map((r.x, r.y)) do d
        (d[1:(end-1)] .+ d[2:end]) ./ 2
    end
    DD.X(xmid), DD.Y(ymid),DD.Dim{:ZONE}(1:60), DD.Dim{:Hemisphere}(1:2)
end

function Base.show(io::IO, tree::UTMTree)
    # Check if this is a compact display (when used within other show methods)
    compact = get(io, :compact, false)

    if compact
        # Compact format for use within other show methods
        print(io, "UTMTree($(length(tree.x))×$(length(tree.y)))x60x2")
    else
        # Full format with copy-pastable constructor on first line
        print(io, "UTMTree($(length(tree.x))×$(length(tree.y)) array, chunksize)")

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
    UTMTree(x, y, transform=UnitSphereFromUTM())

Constructs a UTM Tree with x and y denoting the bounds in each UTM zone.
"""
UTMTree(x, y) = UTMTree(x, y, UnitSphereFromUTM())

"""
    UTMTree(ar::DD.AbstractDimArray,spatial_dims;transform=UnitSphereFromUTM())

Convenience constructor to create a unitspherical spatial search tree from an AbstractDimArray. The spatial dimension
names can be passed as a tuple of symbols. Tries to guess cell boundaries. 
"""
function UTMTree(ar::DD.AbstractDimArray, spatial_dims=(DD.XDim, DD.YDim, DD.Dim{:ZONE}(), DD.Dim{:Hemisphere}()); transform=UnitSphereFromUTM())
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
    xr, yr = map(x->boundrangefromcenters(x.val), ar_spatial_dims[1:2])
    return UTMTree(xr, yr, transform)
end

get_tag(::UTMTree) = nothing
nlevel(r::UTMTree) = max(ceil(Int, log2(length(r.x)))+6, ceil(Int, log2(length(r.y)))+1)
rootnode(t::UTMTree) = TreeNode(t, UTMIndex((1, length(t.x)), (1, length(t.y)),(1,61),3))
extent(t::UTMTree, index::UTMIndex) = Extent(
    X=(t.x[index.x[1]], t.x[index.x[2]]), 
    Y=(t.y[index.y[1]], t.y[index.y[2]]),
    ZONE=(index.zone[1],index.zone[2]-1),
    Hemisphere= index.hemi==3 ? (1,2) : (index.hemi,index.hemi)
)
function linind(grid::UTMTree, index::UTMIndex)
    LinearIndices((length(grid.x)-1, length(grid.y)-1,60,2))[index.x[1], index.y[1],index.zone[1],index.hemi]
end
isleaf(index::UTMIndex) = _isone(index.x) && _isone(index.y) && _isone(index.zone) && (index.hemi !== 3)
function circle_from_extent_1(ex, grid::UTMTree)
    trans=grid.trans
    (x1, x2), (y1, y2), (zone1,zone2),(hemi1,hemi2) = bounds(ex)
    if hemi1 != hemi2
        return SphericalCap(UnitSphereFromGeographic()((0.0,0.0)),π)
    end
    cx, cy, cz = (x2 + x1) / 2, (y2 + y1) / 2, (zone1 + zone2) ÷ 2
    a, b, c, d, e, f, g, h = map(trans, 
        ((x1, y1, zone1, hemi1), (x2, y1, zone2, hemi1), (x2, y2,zone2,hemi1), (x1, y2,zone1, hemi1), 
        (cx, y1, cz, hemi1), (x2, cy, zone2, hemi1), (cx, y2, cz, hemi1), (x1, cy, zone1, hemi1))
    )
    z = trans((cx, cy, cz, hemi1))
    alld = map(p->spherical_distance(z, p), (a, b, c, d, e, f, g, h))
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

function node_to_polygon_unitsphere(grid::UTMTree, index::UTMIndex)
    x1, x2 = index.x
    y1, y2 = index.y
    xr = grid.x
    yr = grid.y
    zone1,zone2 = index.zone
    if index.hemi == 3
        error("Can not make polygon across hemispheres")
    end
    poly = #= @SVector =#[(xr[x1], yr[y1],zone1), (xr[x2], yr[y1], zone2-1), (xr[x2], yr[y2],zone2-1), (xr[x1], yr[y2], zone1), (xr[x1], yr[y1], zone1)]
    grid.trans.(poly)
end
index_to_cartesian(i::Integer, t::UTMTree) = CartesianIndices((length(t.x) - 1, length(t.y) - 1, 60, 2))[i].I

function index_to_native_coords(i, t::UTMTree)
    xhalfstep = get_step(t.x) / 2
    yhalfstep = get_step(t.y) / 2
    i, j, z, h = index_to_cartesian(i, t)
    x = t.x[i] + xhalfstep
    y = t.y[j] + yhalfstep
    x, y, z, h
end

function TreeNode(tree::UTMTree, targetinds::Tuple)
    r1, r2, r3, r4 = targetinds
    ix1, ix2 = first(r1), last(r1)
    iy1, iy2 = first(r2), last(r2)
    TreeNode(tree, UTMIndex((ix1, ix2+1), (iy1, iy2+1), (first(r3),last(r3)+1), (length(r4) == 2 ? 3 : first(r4))))
end

function ProjectionSource(::Type{<:UTMTree}, ar, spatial_dims=(DD.XDim, DD.YDim, DD.Dim{:ZONE},DD.Dim{:Hemisphere}))
    tree = UTMTree(ar, spatial_dims)
    lookups = map(DD.format, DD.dims(ar, spatial_dims))
    chunks = map(eachchunk(ar.data).chunks, DD.dims(ar)) do c, d
        DD.rebuild(d, c)
    end
    xchunks, ychunks = DD.dims(chunks, spatial_dims)
    xchunkbnds = vcat(tree.x[first.(xchunks.val)], last(tree.x))
    ychunkbnds = vcat(tree.y[first.(ychunks.val)], last(tree.y))
    chunktree = UTMTree(xchunkbnds, ychunkbnds)
    ProjectionSource(ar, tree, chunktree, lookups, chunks)
end

#Compute indices given a chunk index for the high-resolution tree
function indices_from_chunk(s::ProjectionSource{<:Any,<:UTMTree}, target_chunk)
    inds = index_to_cartesian(target_chunk, s.chunktree)
    chunkrange = map(getindex, s.chunks, inds)
    map(chunkrange) do cr
        Colon()(extrema(cr)...)
    end
end

function ProjectionTarget(::Type{<:UTMTree}, x, y, trans=UnitSphereFromUTM(); chunksize=512)
    tree = UTMTree(x, y, trans)
    xchunk = x[1:chunksize:end]
    ychunk = y[1:chunksize:end]
    chunktree = UTMTree(xchunk, ychunk, trans)
    ProjectionTarget(tree, chunktree)
end
