import DiskArrays: eachchunk
import DimensionalData as DD
import GeometryOps.SpatialTreeInterface: dual_depth_first_search
using DiskArrays: AbstractDiskArray, findchunk, DiskArrays, ChunkIndex
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

function Base.show(io::IO, ::MIME"text/plain", ps::ProjectionSource)
    T = eltype(ps.ar)
    dims = size(ps.ar)
    dims_str = join(dims, "×")
    print(io, "ProjectionSource{$T}($dims_str, ")
    show(io, ps.tree)
    print(io, ")")
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

function compute_connected_chunks(source::ProjectionSource,target::ProjectionTarget)
    
    connected_chunks = [Int[] for _ in 1:nleaf(target.chunktree)]

    dual_depth_first_search(_intersects, rootnode(target.chunktree), rootnode(source.chunktree)) do n1, n2
        push!(connected_chunks[n1], n2)
    end
    connected_chunks
end

function compute_connected_chunks(source::ProjectionSource, target::ProjectionTarget, targetinds)

    target_smalltree = get_subtree(target.tree,targetinds)
    circle = get_gridextent(target.tree, targetinds...)
    pred = Base.Fix1(_intersects, circle)
    res = Int[]
    depth_first_search(pred, rootnode(source.chunktree)) do n
        test_intersect_highres(source,target_smalltree, n) && push!(res, n)
    end
    res
end

function test_intersect_highres(source,target_smalltree,sourcechunk)
    ssmallinds = indices_from_chunk(source, sourcechunk)
    source_smalltree = get_subtree(source.tree,ssmallinds)
    any_intersect(target_smalltree, source_smalltree)
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

function Base.show(io::IO, ::MIME"text/plain", lpda::LazyProjectedDiskArray{T}) where T
    dims = size(lpda)
    dims_str = join(dims, "×")
    print(io, "$dims_str LazyProjectedDiskArray{$T}")
end

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

struct NearestProjection end

function compute_indices(a::LazyProjectedDiskArray, targetinds, index_arraybuffer)
    source = a.source
    target = a.target
    targettree = target.tree
    isourcetrans = Base.inv(get_projection(source.tree))
    lookups = DD.dims(source.lookups,source.chunks)
    chunks = source.chunks
    compute_nearest_per_chunk(targetinds, targettree, isourcetrans, lookups, chunks, index_arraybuffer)
end







function DiskArrays.readblock!(a::LazyProjectedDiskArray, aout, targetinds::AbstractUnitRange...; index_arraybuffer=make_indexbuffer(a.source.tree, a.target.tree))
    outarray = OffsetArray(aout, targetinds...)
    chunks = compute_connected_chunks(a.source, a.target,targetinds)
    isourcetrans = inv(get_projection(a.source.tree))
    if length(chunks) < 8
        project_batched(a,outarray,chunks,isourcetrans,targetinds)
    else
        project_sequential(a,outarray,chunks,isourcetrans,targetinds;index_arraybuffer)
    end
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

include("sequential.jl")
include("batched.jl")