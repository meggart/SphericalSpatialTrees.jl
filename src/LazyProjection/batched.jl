function project_kernel_batched!(::NearestProjection, outar, targetinds,sourcearrays, targettree, isourcetrans, lookups,chunks)
    alllinind = LinearIndices(gridsize(targettree))
    Threads.@threads for targetindex in CartesianIndices(targetinds)
    #for targetindex in CartesianIndices(targetinds)
        ind = alllinind[targetindex]
        unit = index_to_unitsphere(ind, targettree)
        sourcecoords = isourcetrans(unit)
        sourceindices = map(sourcecoords,lookups) do coord,look
            DD.selectindices(look, DD.Near(coord))
        end
        chunkindices = map((c,i)->findchunk(c.val,i),chunks,sourceindices)
        cI = CartesianIndex(chunkindices)
        outar[targetindex] = sourcearrays[cI][sourceindices...]
    end
    outar
end

function load_sourcechunks(source,chunks)
    input_chunks = map(chunks) do c
        chunk_cartindex = CartesianIndex(index_to_cartesian(c,source.chunktree))
        sourceindices = indices_from_chunk(source,c)
        aout = source.ar.data[sourceindices...]
        data = OffsetArray(aout,sourceindices...)
        chunk_cartindex => data
    end
    nsource = ndims(source.ar.data)
    i1,i2 = extrema(first,input_chunks)
    s = (i2.-i1+oneunit(i1)).I
    et = typeof(last(first(input_chunks)))
    sourcearrays = OffsetArray(Array{et,nsource}(undef,s...),(i1-oneunit(i1)).I...)
    for (i,a) in input_chunks
        sourcearrays[i]=a
    end
    sourcearrays
end

function project_batched(a::LazyProjectedDiskArray,outarray,chunks,isourcetrans,targetinds)
    sourcearrays = load_sourcechunks(a.source,chunks)
    project_kernel_batched!(NearestProjection(),outarray,targetinds,sourcearrays, a.target.tree, isourcetrans, a.source.lookups,a.source.chunks)
end