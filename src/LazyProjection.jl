import DiskArrays: eachchunk
import DimensionalData as DD
import GeometryOps.SpatialTreeInterface: do_dual_query
using DiskArrays: AbstractDiskArray, findchunk, DiskArrays
using OffsetArrays: OffsetArray
using ProgressMeter
using Distributed: pmap

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

function compute_connected_chunks(source::ProjectionSource, target::ProjectionTarget, targetinds)

    circle = get_gridextent(target.tree, targetinds...)
    pred = Base.Fix1(_intersects, circle)
    res = Int[]
    do_query(pred, rootnode(source.chunktree)) do n
        push!(res, n)
    end
    res
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

function compute_nearest_per_chunk(targetinds, targettree, isourcetrans, lookups::Tuple{Vararg{<:Any,Nsource}}, chunks, index_arraybuffer) where Nsource
    alllinind = LinearIndices(gridsize(targettree))
    #Ntarget = ndims(targettree)
    inner_indexarray = fill((zero(CartesianIndex{Nsource}), zero(CartesianIndex{Nsource})), length.(targetinds)...)
    indexarray = OffsetArray(inner_indexarray, targetinds...)
    Threads.@threads for targetindex in CartesianIndices(targetinds)
        ind = alllinind[targetindex]
        unit = index_to_unitsphere(ind, targettree)
        sourcecoords = isourcetrans(unit)
        sourceindices = map(sourcecoords,lookups) do coord,look
            DD.selectindices(look, DD.Near(coord))
        end
        chunkindices = map((c,i)->findchunk(c.val,i),chunks,sourceindices)
        cI = CartesianIndex(chunkindices)
        indexarray[targetindex] = (cI, CartesianIndex(sourceindices))
    end
    cartinds = first.(unique(first, indexarray))
    if length(cartinds) > length(index_arraybuffer)
        error("Too many connected chunks")
    end
    mybuffer = view(index_arraybuffer, 1:length(cartinds))
    foreach(mybuffer) do b
        empty!(first(b))
        empty!(last(b))
    end
    for itarget in CartesianIndices(indexarray)
        chunknum, iel = indexarray[itarget]
        ichunk = findfirst(==(chunknum), cartinds)
        vt, vs = index_arraybuffer[ichunk]
        push!(vt, itarget)
        push!(vs, iel)
    end
    return mybuffer
end

function compute_indices(a::LazyProjectedDiskArray, targetinds, index_arraybuffer)
    source = a.source
    target = a.target
    targettree = target.tree
    isourcetrans = Base.inv(get_projection(source.tree))
    lookups = DD.dims(source.lookups,source.chunks)
    chunks = source.chunks
    compute_nearest_per_chunk(targetinds, targettree, isourcetrans, lookups, chunks, index_arraybuffer)
end

function copydata(outarray, inds, source)
    for (vt, vs) in inds
        i1,i2 = extrema(vs)
        bbr = map(Colon(),i1.I,i2.I)
        data = OffsetArray(source[bbr...], bbr...)
        outarray[vt] = data[vs]
    end
end

function make_indexbuffer(sourcetree, targettree, N=50)
    Nsource = ndims(sourcetree)
    Ntarget = ndims(targettree)
    [(CartesianIndex{Ntarget}[], CartesianIndex{Nsource}[]) for _ in 1:N]
end

function DiskArrays.readblock!(a::LazyProjectedDiskArray, aout, targetinds::AbstractUnitRange...; index_arraybuffer=make_indexbuffer(a.source.tree, a.target.tree))
    outarray = OffsetArray(aout, targetinds...)
    inds = compute_indices(a, targetinds, index_arraybuffer)
    copydata(outarray, inds, a.source.ar.data)
end

function reproject!(target_array,source,target)
    #this assumes there are only x and y axis
    lazyarray = LazyProjectedDiskArray(source,target)
    targetchunks = eachchunk(target_array)
    index_arraybuffer = make_indexbuffer(source.tree, target.tree)
    aout = zeros(eltype(target_array.data), length.(first(targetchunks))...)
    @showprogress for targetchunk in targetchunks
        DiskArrays.readblock!(lazyarray, aout, targetchunk...; index_arraybuffer)
        target_array.data[targetchunk...] = aout
    end
end