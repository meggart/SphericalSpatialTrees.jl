import DiskArrays: eachchunk
import DimensionalData as DD
import GeometryOps.SpatialTreeInterface: do_dual_query
using DiskArrays: AbstractDiskArray, findchunk, DiskArrays
using OffsetArrays: OffsetArray

#Some helper functions to detect the bounds of an area from their centers. 
# This needs much better integration with DD
function mayberound(x,st)
    xrnd = round(x,digits=-floor(Int,log10(abs(st))))
    (xrnd-x)/x < 1e-9 ? xrnd : x
end
function boundrangefromcenters(x)
    lb = mayberound(first(x) - step(x)/2, step(x))
    ub = mayberound(last(x)+step(x)/2,step(x))
    range(lb,ub,length(x)+1)
end


struct ProjectionSource{Y<:DD.AbstractDimArray,T,L,C,CT}
    ar::Y
    tree::T
    chunktree::CT
    lookups::L
    chunks::C
end
function ProjectionSource(::Type{<:RegularGridTree}, ar,spatial_dims = (:X,:Y))
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

function compute_connected_chunks(source::ProjectionSource,target::ProjectionTarget)
    
    connected_chunks = [Int[] for _ in 1:nleaf(target.chunktree)]

    do_dual_query(_intersects, rootnode(target.chunktree), rootnode(source.chunktree)) do n1, n2
        push!(connected_chunks[n1], n2)
    end
    connected_chunks
end

struct LazyProjectedDiskArray{T,N,S,TA} <: AbstractDiskArray{T,N}
    source::ProjectionSource
    target::ProjectionTarget
end
function LazyProjectedDiskArray(source,target)
    LazyProjectedDiskArray{eltype(source.ar),ndims(target.tree),typeof(source),typeof(target)}(source,target)
end
Base.size(a::LazyProjectedDiskArray) = gridsize(a.target.tree)
Base.ndims(a::LazyProjectedDiskArray) = ndims(a.target.tree)

function compute_nearest_per_chunk(targetinds, targettree, isourcetrans,lookups,chunks)
    alllinind = LinearIndices(gridsize(targettree))
    Ntarget = ndims(targettree)
    Nsource = length(lookups)
    sourced = Dict{CartesianIndex{Nsource},Tuple{Vector{CartesianIndex{Nsource}},Vector{CartesianIndex{Ntarget}}}}()
    for targetindex in CartesianIndices(targetinds)
        ind = alllinind[targetindex]
        unit = index_to_unitsphere(ind, targettree)
        sourcecoords = isourcetrans(unit)
        sourceindices = map(sourcecoords,lookups) do coord,look
            DD.selectindices(look, DD.Near(coord))
        end
        chunkindices = map((c,i)->findchunk(c.val,i),chunks,sourceindices)
        cI = CartesianIndex(chunkindices)
        vs,vt = get!(()->(CartesianIndex{Nsource}[],CartesianIndex{Ntarget}[]),sourced,cI)
        push!(vt,targetindex)
        push!(vs,CartesianIndex(sourceindices))
    end
    sourced
end


function DiskArrays.readblock!(a::LazyProjectedDiskArray,aout,targetinds::AbstractUnitRange...)
    source = a.source
    target = a.target
    targettree = target.tree
    isourcetrans = Base.inv(get_projection(source.tree))
    lookups = DD.dims(source.lookups,source.chunks)
    chunks = source.chunks
    inds = compute_nearest_per_chunk(targetinds, targettree, isourcetrans,lookups,chunks)
    outarray = OffsetArray(aout,targetinds...) 
    for (_,(vs,vt)) in inds
        i1,i2 = extrema(vs)
        bbr = map(Colon(),i1.I,i2.I)
        data = OffsetArray(source.ar.data[bbr...],bbr...)
        outarray[vt] = data[vs]
    end
end