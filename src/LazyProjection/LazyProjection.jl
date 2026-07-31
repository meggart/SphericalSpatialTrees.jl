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
    st = pseudostep(x)
    lb = mayberound(first(x) - st / 2, st)
    ub = mayberound(last(x) + st / 2, st)
    range(lb,ub,length(x)+1)
end
pseudostep(a::AbstractRange) = step(a)
pseudostep(a::AbstractVector) = length(a) > 1 ? (last(a) - first(a)) / (length(a)-1) : one(eltype(a))


"""
    ProjectionSource(::Type{<:SpatialTree}, ar, spatial_dims...)

A regridding source: a chunked array `ar` together with a spatial search tree
(`tree`), a coarser chunk tree (`chunktree`) whose leaves correspond to the
chunks of `ar`, and the lookup axes (`lookups`) used for nearest-neighbor
searches.

Construct it by passing the tree type as the first argument:

    source = ProjectionSource(RegularGridTree, geo_array)

The constructor is defined per tree type; see the documentation of the
respective tree for its specific arguments (`spatial_dims`, etc.).
"""
struct ProjectionSource{Y<:DD.AbstractDimArray,T,L,C,CT}
    ar::Y
    tree::T
    chunktree::CT
    lookups::L
    chunks::C
end

function Base.show(io::IO, ::MIME"text/plain", ps::ProjectionSource)
    # First line: Show tree constructor in copy-pastable format
    b = IOBuffer()
    bcompact = IOContext(b, :compact => true)
    show(bcompact,ps.tree)
    treestring = String(take!(b))
    
    # Additional lines in cyan: Show ProjectionSource{T}(dims) with element type and dimensions
    T = eltype(ps.ar)
    dims = size(ps.ar)
    dims_str = join(dims, "×")
    printstyled(io, "ProjectionSource{$T}($dims_str, $treestring)", color=:cyan)
end

"""
    ProjectionTarget(::Type{<:SpatialTree}, args...; kwargs...)

A regridding target: a spatial search tree (`tree`) together with a coarser
chunk tree (`chunktree`) that defines the chunks of the resulting array.

Construct it by passing the tree type as the first argument:

    target = ProjectionTarget(RegularGridTree, -180.0:0.1:180.0, 90.0:-0.1:-90.0)
    target = ProjectionTarget(ISEACircleTree, 8, 2)

The constructor is defined per tree type; see the documentation of the
respective tree for its specific arguments (resolutions, `chunksize`, etc.).
"""
struct ProjectionTarget{T,CT}
    tree::T
    chunktree::CT
end
function Base.show(io::IO, ::MIME"text/plain", ps::ProjectionTarget)
    # First line: Show tree constructor in copy-pastable format
    b = IOBuffer()
    bcompact = IOContext(b, :compact => true)
    show(bcompact,ps.tree)
    treestring = String(take!(b))

    printstyled(io, "ProjectionTarget($treestring)", color=:cyan)
end

"""
    create_dataset(target, path; arrayname=:layer, arraymeta=Dict(), datasetmeta=Dict(),
                   backend=:zarr, output_datatype=Float64, kwargs...)

Create a dataset on disk at `path` whose grid matches the target tree of
`target`. It contains a single array `arrayname` that is chunked according to
the target's chunk tree. Returns the created array, ready to be filled by
[`reproject!`](@ref).
"""
function create_dataset(target::ProjectionTarget, 
    path; arrayname=:layer, arraymeta=Dict(), datasetmeta=Dict(), backend=:zarr, output_datatype=Float64, kwargs...)
    
    back = YAXArrayBase.backendfrompath(path;driver=backend)
    dims = DD.dims(target.tree)
    chunkdims = DD.dims(target.chunktree)
    cs = map(dims,chunkdims) do d,cd
        length(d) ÷ length(cd)
    end
    dimnames = string.(DD.name.(chunkdims))
    group = YAXArrayBase.create_dataset(
            back,
            path,
            datasetmeta,
            dimnames,
            map(d->d.val,dims),
            (Dict(),Dict(),Dict()),
            (output_datatype,),
            (string(arrayname),),
            (dimnames,),
            (arraymeta,),
            (cs,);
            kwargs...
        )
    group[string(arrayname)]
end


"""
    compute_connected_chunks(source, target)
    compute_connected_chunks(source, target, targetinds)

Determine which source chunks are needed for regridding. The two-argument
version returns a vector with one entry per target chunk (in linear index
order), each entry holding the indices of the source chunks that intersect it.
The three-argument version restricts the computation to the given range of
target indices and returns only the source chunk indices needed there.
"""
function compute_connected_chunks(source::ProjectionSource,target::ProjectionTarget)
    
    connected_chunks = [Int[] for _ in 1:nleaf(target.chunktree)]

    dual_depth_first_search(_intersects, rootnode(target.chunktree), rootnode(source.chunktree)) do n1, n2
        push!(connected_chunks[n1], n2)
    end
    connected_chunks
end

function compute_connected_chunks(source::ProjectionSource, target::ProjectionTarget, targetinds)
    res = Int[]
    with_transform(target.tree) do targettree
        with_transform(source.chunktree) do sourcechunktree
            with_transform(source.tree) do sourcetree
                target_smalltree = TreeNode(targettree, targetinds)
                circle = get_gridextent(targettree, targetinds...)
                pred = Base.Fix1(_intersects, circle)
                depth_first_search(pred, rootnode(sourcechunktree)) do n
                    test_intersect_highres(source, target_smalltree, n, sourcetree) && push!(res, n)
                end
            end
        end
    end
    res
end

function test_intersect_highres(source, target_smalltree, sourcechunk, sourcetree)
    ssmallinds = indices_from_chunk(source, sourcechunk)
    source_smalltree = TreeNode(sourcetree, ssmallinds)
    any_intersect(target_smalltree, source_smalltree)
end


"""
    LazyProjectedDiskArray(source, target)

A lazy `AbstractDiskArray` that regrids data from `source` to `target` by
nearest-neighbor search, computed on demand. It has the size of the target grid
and is chunked according to the target's chunk tree; accessing a block triggers
loading and reprojection of the connected source chunks.
"""
struct LazyProjectedDiskArray{T,N,S,TA} <: AbstractDiskArray{T,N}
    source::ProjectionSource
    target::ProjectionTarget
end
function LazyProjectedDiskArray(source,target)
    LazyProjectedDiskArray{eltype(source.ar),ndims(target.tree),typeof(source),typeof(target)}(source,target)
end
function DiskArrays.eachchunk(a::LazyProjectedDiskArray)
    gs = gridsize(a.target.tree)
    cgs = gridsize(a.target.chunktree)
    cs = Int.(gs .÷ cgs)
    DiskArrays.GridChunks(a,cs)
end
DiskArrays.haschunks(::LazyProjectedDiskArray) = DiskArrays.Chunked()
Base.size(a::LazyProjectedDiskArray) = gridsize(a.target.tree)
Base.ndims(a::LazyProjectedDiskArray) = ndims(a.target.tree)

DD.DimArray(source::ProjectionSource, target::ProjectionTarget) = 
    DD.DimArray(LazyProjectedDiskArray(source,target),DD.dims(target.tree))



function Base.show(io::IO, ::MIME"text/plain", lpda::LazyProjectedDiskArray{T}) where T
    dims = size(lpda)
    dims_str = join(dims, "×")
    print(io, "$dims_str LazyProjectedDiskArray{$T}")
end

function compute_nearest_per_chunk(targetinds, targettree, isourcetrans, lookups::Tuple{Vararg{Any,Nsource}}, chunks, index_arraybuffer) where Nsource
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
    if length(chunks) < 100
        project_batched(a,outarray,chunks,isourcetrans,targetinds)
    else
        project_sequential(a,outarray,chunks,isourcetrans,targetinds;index_arraybuffer)
    end
end

"""
    reproject!(target_array, source, target)

Regrid all data from `source` to `target`, writing the result chunk by chunk
into `target_array` (e.g. a Zarr array created with [`create_dataset`](@ref)).
This assumes that `target_array` only has spatial axes.
"""
function reproject!(target_array,source,target)
    #this assumes there are only spatial axes
    lazyarray = LazyProjectedDiskArray(source,target)
    targetchunks = eachchunk(target_array)
    index_arraybuffer = make_indexbuffer(source.tree, target.tree)
    aout = zeros(eltype(target_array), length.(first(targetchunks))...)
    @showprogress for targetchunk in targetchunks
        DiskArrays.readblock!(lazyarray, aout, targetchunk...; index_arraybuffer)
        target_array[targetchunk...] = aout
    end
end

include("threading_helpers.jl")
include("sequential.jl")
include("batched.jl")