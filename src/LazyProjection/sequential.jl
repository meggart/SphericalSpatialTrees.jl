function copy_sequential(outarray, inds, source)
    for (vt, vs) in inds
        i1,i2 = extrema(vs)
        bbr = map(Colon(),i1.I,i2.I)
        data = OffsetArray(source[bbr...], bbr...)
        outarray[vt] = data[vs]
    end
end

function project_sequential(a::LazyProjectedDiskArray,outarray,chunks,isourcetrans,targetinds)
    indices_per_chunk = precompute_sequential_weights(targetinds, a.target.tree, isourcetrans, a.source.lookups, s.source.chunks, index_arraybuffer)
    copy_sequential(outarray,indices_per_chunk,a.source.ar)
end


function precompute_sequential_weights(targetinds, targettree, isourcetrans, lookups::Tuple{Vararg{<:Any,Nsource}}, chunks, index_arraybuffer) where Nsource
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
