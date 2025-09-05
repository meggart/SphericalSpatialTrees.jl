import DiskArrays: eachchunk
import DimensionalData as DD
import GeometryOps.SpatialTreeInterface: do_dual_query

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