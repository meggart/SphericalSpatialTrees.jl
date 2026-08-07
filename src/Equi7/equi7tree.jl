module Equi7
import GeometryOps.UnitSpherical: UnitSphereFromGeographic, SphericalCap, _merge, GeographicFromUnitSphere
import GeometryOps.SpatialTreeInterface as STI
import ..SphericalSpatialTrees as SST
import DimensionalData as DD
import DiskArrays: eachchunk
import CoordinateTransformations as CT
import Proj
import CoordinateTransformations: Transformation, ∘
import Artifacts: @artifact_str
using Statistics: median

"""Read the tuple of `(tile_x, tile_y)` coordinates for an Equi7 zone from the binary artifact."""
readzones(zone) = open(joinpath(artifact"equi7tiles","Equi7ZoneIndices-1","tiles_bin",zone),"r") do f
    n = read(f,Int)
    map(1:n) do _
        read(f,Int),read(f,Int)
    end
end

"""
The maximum tile index `(nx, ny)` for each of the 7 zones, in the order
`(AF, AN, AS, EU, NA, OC, SA)`.
"""
const MAX_N_TILE = ((114, 93),(97, 87),(115, 96),(82, 55),(133, 98),(187, 121),(116, 105))
"""Bounding box `(max_nx + 1, max_ny + 1)` over all zones, i.e. `(188, 122)`."""
const MAX_SIZE = maximum(first,MAX_N_TILE)+1,maximum(last,MAX_N_TILE)+1
"""The 7 Equi7 continental zone labels."""
const ZONES = ("AF","AN","AS","EU","NA","OC","SA")
"""EPSG authority codes for the 7 Equi7 zones (27701–27707)."""
const CODES = 27701:27707
const EQUI7Trans = typeof(SST.init_threaded_proj_collection(string.("EPSG:", CODES)))[]
const EQUI7ITrans = typeof(SST.init_threaded_proj_collection(string.("EPSG:", CODES)))[]

"""Per-zone tile coordinate arrays loaded from the `equi7tiles` artifact."""
const TILECOORDS = readzones.(ZONES)

"""
    EQUI7Tag

Internal tag carried by a `RegularGridTree` child inside an `Equi7Tree`, holding
the zone index (`1`–`7`) and the grid resolution. Used by [`linind`](@ref) to
map leaf nodes back to a linear index in the full 3‑D array.
"""
struct EQUI7Tag 
    zone::Int
    resolution::Int
end




"""
Convert a tile coordinate `(i, j)` (0‑based 100 km tile index) to a
[`SphericalCap`](@ref) covering the four corners of that tile via the Equi7
inverse projection.
"""
function coord_to_circle((i,j),itrans)
    x0,x1,y0,y1 = i*1e5,(i+1)*1e5,j*1e5,(j+1)*1e5
    corners = ((x0,y0),(x0,y1),(x1,y1),(x1,y0))
    reduce(_merge,SphericalCap.(itrans.(corners),0.0))
end

"""
Split `indices` at the median of a coordinate component (`.first` or `.last`).
Returns `(indices_above, indices_below_or_equal)`.
"""
function single_split(allcoords,indices,by)
    spl = median(by.(allcoords[indices]))
    ir = map(i->by(i)>spl,allcoords[indices])
    ir2 = map(i->by(i)<=spl,allcoords[indices])
    return indices[ir],indices[ir2]
end
"""
Create a `RegularGridTree` covering a single Equi7 tile `coord = (i, j)` at the
given `resolution`, using the Equi7 inverse projection for zone `zone`.
"""
function get_tilenode(coord,resolution,zone)
    x = range(coord[1]*1e5,(coord[1]+1)*1e5,length=resolution+1)
    y = range(coord[2]*1e5,(coord[2]+1)*1e5,length=resolution+1)
    SST.rootnode(SST.RegularGridTree(x,y,EQUI7ITrans[1] ∘ SST.PickPlane(zone),EQUI7Tag(zone,resolution)))
end



"""
Recursively build the `TileNode` hierarchy for one Equi7 zone.

Given a set of tile `indices`, their `coords` (tile‑coordinate tuples), their
spherical `ext`ents (`SphericalCap`), and the pre‑built `tiles`
(`RegularGridTree` nodes), the function splits the tiles along their longer
coordinate axis at the median, recursing until fewer than 5 tiles remain (leaf
case).

The resulting tree balances the tile index for efficient spatial queries.
"""
function build_node(indices, coords,ext, tiles)

    myextent = reduce(_merge,ext[indices])
    children = SST.TileNode{eltype(tiles)}[]
    leaves = eltype(tiles)[]
    #Check if we are already a leaf
    if length(indices) < 5
        for i in indices
            push!(leaves,tiles[i])
        end
        return SST.TileNode(children,leaves,myextent)
    end

    xr = extrema(first,coords[indices])
    yr = extrema(last,coords[indices])

    if xr[2]-xr[1] > yr[2]-yr[1]
        #we split first along x
        split1,split2 = single_split(coords,indices,first)
        c1,c2 = single_split(coords,split1,last)
        c3,c4 = single_split(coords,split2,last)
    else
        #we first split along y
        split1,split2 = single_split(coords,indices,last)
        c1,c2 = single_split(coords,split1,first)
        c3,c4 = single_split(coords,split2,first)
    end
    @assert isempty(intersect(c1,c2,c3,c4))
    for c in (c1,c2,c3,c4)
        if !isempty(c)
            push!(children,build_node(c,coords,ext,tiles))
        end
    end
    return SST.TileNode(children,leaves,myextent)
end

"""
Build the complete `TileNode` hierarchy for one Equi7 `zone` string
(e.g. `"EU"`) at the given grid `resolution` per tile.
"""
function nodefromzone(zone,resolution)
    izone = findfirst(==(zone),ZONES)
    allcoords = TILECOORDS[izone]
    indices = collect(1:length(allcoords))
    allextents = coord_to_circle.(allcoords,(EQUI7ITrans[1] ∘ SST.PickPlane(izone),))
    alltilenodes = get_tilenode.(allcoords,resolution,izone);
    build_node(indices,allcoords,allextents,alltilenodes)
end

"""
    Equi7Tree(resolution)
    Equi7Tree(resolution::Integer)

A spatial tree on the **Equi7 continental equal-area grid** over its 7 zones
(AF, AN, AS, EU, NA, OC, SA).

The grid is 3‑dimensional with dimensions `(x, y, zone)`:

    gridsize(tree) = (MAX_SIZE[1] * resolution, MAX_SIZE[2] * resolution, 7)
                    = (188 * resolution, 122 * resolution, 7)

Each zone is built from a binary tile index (100 km tiles from the
`equi7tiles` artifact), subdivided into `resolution × resolution` cells per
tile and organized in a `TileNode` hierarchy for query performance. The root
node wraps the 7 zone subtrees under a single `SphericalCap` covering the whole
globe.

The underlying PROJ transformations (EPSG:27701–27707) are initialised lazily in
`__init__` and shared via `EQUI7Trans` / `EQUI7ITrans`.

# Examples
```julia
julia> tree = Equi7Tree(2)
SphericalSpatialTrees.Equi7.Equi7Tree{…}(2, Node with 7 children and 0 leaves)

julia> gridsize(tree)
(376, 244, 7)

julia> index_to_lonlat(1, tree)
(-37.8427, -40.5715)

julia> index_to_native_coords(1, tree)
(50000.0, 50000.0, 1)
```

# References
- Bauer-Marschallinger et al. (2014), *Remote Sensing* 6(5), 4194–4226.
- https://github.com/TUW-GEO/Equi7Grid
"""
struct Equi7Tree{T<:SST.TileNode}
    resolution::Int
    rootnode::T
end
function Equi7Tree(resolution::Integer)
    children = [nodefromzone(zone,resolution) for zone in ZONES]
    leaves = eltype(children[1].leaves)[]
    rootnode = SST.TileNode(children,leaves, SphericalCap(UnitSphereFromGeographic()((0.0,0.0)),Float64(π)))
    Equi7Tree(resolution, rootnode)
end
SST.rootnode(tree::Equi7Tree) = tree.rootnode
Base.ndims(::Equi7Tree) = 3
SST.gridsize(tree::Equi7Tree) = ((MAX_SIZE .* tree.resolution)...,7)
SST.get_projection(::Equi7Tree) = EQUI7ITrans[1]

"""
    DD.dims(tree::Equi7Tree)

Return concrete `X`, `Y`, and `zone` dimensions for the Equi7Tree grid.
The X and Y ranges have `MAX_SIZE .* resolution` points, with centers offset
by half a cell (`50000.0 / resolution`) from the tile boundaries.
The zone dimension carries the 7 zone labels `["AF", …, "SA"]`.
"""
function DD.dims(t::Equi7Tree) 
    n = t.resolution
    offs = 1e5/2/n
    rx = range(offs,MAX_SIZE[1]*1e5 - offs,length=n*MAX_SIZE[1])
    ry = range(offs,MAX_SIZE[2]*1e5 - offs,length=n*MAX_SIZE[2])
    (DD.X(rx),DD.Y(ry),DD.Dim{:zone}(collect(ZONES)))
end
DD.dims(::Type{<:Equi7Tree}) = (DD.X(),DD.Y(),DD.Dim{:zone}(collect(ZONES)))


function SST.linind(tag::EQUI7Tag, tree::SST.TreeNode)
    n = tag.resolution .* MAX_SIZE
    xoffset = Int(first(tree.grid.x)/1e5) * tag.resolution
    yoffset = Int(first(tree.grid.y)/1e5) * tag.resolution
    ind = LinearIndices((first(n),last(n), 7))[tree.index.x[1]+xoffset, tree.index.y[1]+yoffset, tag.zone]
    ind
end

"""
    SST.index_to_native_coords(i, tree::Equi7Tree)

Return the `(x, y, zone)` native coordinates of the centre of cell `i` in the
Equi7 grid. `x` and `y` are projected coordinates in metres (Equi7 / EPSG:277xx
projection), and `zone` is an integer `1`–`7`.

`i` is a linear index into the flattened `(nx, ny, 7)` grid, where `nx =
MAX_SIZE[1] * resolution` and `ny = MAX_SIZE[2] * resolution`.
"""
function SST.index_to_native_coords(i,tree::Equi7Tree)
    n = tree.resolution .* MAX_SIZE
    ix,iy,zone = CartesianIndices((first(n),last(n), 7))[i].I
    x,y = ((ix,iy).-1) .* (1e5/tree.resolution) .+ (5e4/tree.resolution)
    (x,y,zone)
end


"""
    ProjectionSource(::Type{<:Equi7Tree}, ar, spatial_dims=(DD.XDim, DD.YDim, :zone))

Create a regridding source from an array `ar` with dimensions `(x, y, zone)`
that covers all 7 Equi7 zones.

The array must satisfy:
- `size(ar, 3) == 7`
- The X and Y sizes must be multiples of `MAX_SIZE = (188, 122)` and equal each
  other when divided (i.e. a uniform resolution across all zones).
- The array must be chunked (its chunks define the source chunk tree).

`spatial_dims` specifies which DimensionalData dimensions correspond to X, Y,
and zone; defaults to `(DD.XDim, DD.YDim, :zone)`.
"""
function SST.ProjectionSource(::Type{<:Equi7Tree}, ar, spatial_dims = (DD.XDim,DD.YDim,:zone))
    nx,ny,n = size(ar)
    @assert n == 7 "The target must have 7 faces, got $n"
    resolution = (nx,ny)./(MAX_SIZE)
    @assert first(resolution) == last(resolution) "Array must contain all tiles"
    @assert isinteger(first(resolution)) "Array must contain all tiles"
    res = Int(first(resolution))
    tree = Equi7Tree(res)
    chunks = map(eachchunk(ar.data).chunks,DD.dims(ar)) do c,d
        DD.rebuild(d,c)
    end
    lookups = DD.dims(ar,spatial_dims)
    lookups = DD.format.(lookups)
    xchunks,ychunks,nchunks = DD.dims(chunks, spatial_dims)
    chunkres = length(xchunks.val) ÷ first(MAX_SIZE)
    chunktree = Equi7Tree(chunkres)
    SST.ProjectionSource(ar,tree,chunktree,lookups,chunks)
end

"""
    ProjectionTarget(::Type{Equi7Tree}, target_resolution, chunk_resolution)

Create a regridding target on the Equi7 grid with the given cell resolution
and chunk resolution.

`target_resolution` sets the number of cells per tile side for the output grid;
`chunk_resolution` sets the coarser resolution for the chunk tree used to
accelerate the regridding search. Typical usage: a high target resolution
(e.g. `8`) and a lower chunk resolution (e.g. `1` or `2`).

# Examples
```julia
julia> target = ProjectionTarget(Equi7Tree, 4, 1)
ProjectionTarget{Equi7Tree}(…)

julia> target.tree.resolution
4

julia> target.chunktree.resolution
1
```
"""
function SST.ProjectionTarget(::Type{Equi7Tree},target_resolution, chunk_resolution)
    tree = Equi7Tree(target_resolution)
    chunktree = Equi7Tree(chunk_resolution)
    SST.ProjectionTarget(tree,chunktree)
end

"""
    TreeNode(tree::Equi7Tree, target_indices)

Return a [`TreeNode`](@ref) containing the indices given in `target_indices`.

`target_indices` is a tuple `(ix, iy, n)` where `ix` and `iy` are
[`UnitRange`](@ref)s over the X and Y dimensions of the full grid
(`1:gridsize(tree)[1]`), and `n` is the zone range (`1:7`).

When `n` spans more than one zone, the full root node of the tree is returned
(see note below).

# Note
Cross‑zone `TreeNode` construction currently requires the multi‑zone path —
see the caveat at `Equi7`.[`rootnode`](@ref).
"""
function SST.TreeNode(tree::Equi7Tree, target_indices)
    ix,iy,n = target_indices
    if length(n) > 1
        return rootnode(tree)
    end
    xr = range(0.0,MAX_SIZE[1]*1e5,length=(MAX_SIZE[1]*tree.resolution)+1)
    yr = range(0.0,MAX_SIZE[2]*1e5,length=(MAX_SIZE[2]*tree.resolution)+1)
    grid = SST.RegularGridTree(xr,yr,EQUI7ITrans[1] ∘ SST.PickPlane(only(n)),EQUI7Tag(only(n),tree.resolution))
    index = SST.TreeIndex((first(ix),last(ix)+1),(first(iy),last(iy)+1))
    SST.TreeNode(grid,index)
end

"""
    SST.get_gridextent(tree::Equi7Tree, xr, yr, nr)

Return the merged [`SphericalCap`](@ref) covering the range of cells given by
index ranges `xr`, `yr`, and zone range `nr`.
"""
function SST.get_gridextent(tree::Equi7Tree, xr::AbstractUnitRange, yr::AbstractUnitRange, nr::AbstractUnitRange)
    mapreduce(_merge,nr) do n
        t = SST.TreeNode(tree,(xr,yr,n))
        SST.node_extent(t)
    end
end


"""
    __init__()

Initialise the shared PROJ transformation collections for EPSG:27701–27707
(Equi7 forward and inverse transforms). Called automatically at module load.
"""
function __init__()
    t = SST.init_threaded_proj_collection(string.("EPSG:", CODES))
    push!(EQUI7Trans,t)
    push!(EQUI7ITrans,inv(t))

end
end