module Equi7
import GeometryOps.UnitSpherical: UnitSphereFromGeographic, SphericalCap, _merge, GeographicFromUnitSphere
import GeometryOps.SpatialTreeInterface as STI
import ..SphericalSpatialTrees as SST
import DimensionalData as DD
import CoordinateTransformations as CT
import Proj
import CoordinateTransformations: Transformation, ∘
using Statistics: median

const MAX_N_TILE = ((114, 93),(97, 87),(115, 96),(82, 55),(133, 98),(187, 121),(116, 105))
const MAX_SIZE = maximum(first,MAX_N_TILE)+1,maximum(last,MAX_N_TILE)+1
const ZONES = ("AF","AN","AS","EU","NA","OC","SA")
const CODES = 27701:27707
const EQUI7Trans = typeof(SST.init_threaded_proj_collection_epsg(CODES))[]
const EQUI7ITrans = typeof(SST.init_threaded_proj_collection_epsg(CODES))[]
const TILECOORDS = [include("tiles/$zone") for zone in ZONES]
struct EQUI7Tag 
    zone::Int
    resolution::Int
end


function __init__()
    t = SST.init_threaded_proj_collection_epsg(CODES)
    push!(EQUI7Trans,t)
    push!(EQUI7ITrans,inv(t))
end

function coord_to_circle((i,j),itrans)
    x0,x1,y0,y1 = i*1e5,(i+1)*1e5,j*1e5,(j+1)*1e5
    corners = ((x0,y0),(x0,y1),(x1,y1),(x1,y0))
    reduce(_merge,SphericalCap.(itrans.(corners),0.0))
end

function single_split(allcoords,indices,by)
    spl = median(by.(allcoords[indices]))
    ir = map(i->by(i)>spl,allcoords[indices])
    ir2 = map(i->by(i)<=spl,allcoords[indices])
    return indices[ir],indices[ir2]
end
function get_tilenode(coord,resolution,zone)
    x = range(coord[1]*1e5,(coord[1]+1)*1e5,length=resolution+1)
    y = range(coord[2]*1e5,(coord[2]+1)*1e5,length=resolution+1)
    SST.rootnode(SST.RegularGridTree(x,y,EQUI7ITrans ∘ SST.PickPlane(zone),EQUI7Tag(zone,resolution)))
end



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

function nodefromzone(zone,resolution)
    izone = findfirst(==(zone),ZONES)
    allcoords = TILECOORDS[izone]
    indices = collect(1:length(allcoords))
    allextents = coord_to_circle.(allcoords,(EQUI7ITrans[1] ∘ SST.PickPlane(izone),))
    alltilenodes = get_tilenode.(allcoords,resolution,izone);
    build_node(indices,allcoords,allextents,alltilenodes)
end

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

function SST.index_to_native_coords(i,tree::Equi7Tree)
    n = tree.resolution .* MAX_SIZE
    ix,iy,zone = CartesianIndices((first(n),last(n), 7))[i].I
    x,y = ((ix,iy).-1) .* (1e5/tree.resolution) .+ (5e4/tree.resolution)
    (x,y,zone)
end


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

function SST.ProjectionTarget(::Type{Equi7Tree},target_resolution, chunk_resolution)
    tree = Equi7Tree(target_resolution)
    chunktree = Equi7Tree(chunk_resolution)
    SST.ProjectionTarget(tree,chunktree)
end

"""
    TreeNode(tree, target_indices)

Returns a Tree node that contains the indices given in target_indices.
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

function SST.get_gridextent(tree::Equi7Tree, xr::AbstractUnitRange, yr::AbstractUnitRange, nr::AbstractUnitRange)
    mapreduce(_merge,nr) do n
        t = SST.TreeNode(tree,(xr,yr,n))
        SST.node_extent(t)
    end
end



end