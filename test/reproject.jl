using Test
import SphericalSpatialTrees as SST
import DimensionalData as DD
import DiskArrays
import Zarr

function make_testarray(::Type{<:SST.ISEACircleTree})
    a = DiskArrays.mockchunks(rand(4,4,10),(2,2,1))
    d = DD.dims(SST.ISEACircleTree(2))
    DD.DimArray(a,d)
end
function make_testarray(::Type{<:SST.RegularGridTree})
    a = DiskArrays.mockchunks(rand(12,6),(5,5))
    d = DD.Dim{:lon}(-165.0:30.0:165.0),DD.Dim{:lat}(75.0:-30.0:-75.0)
    DD.DimArray(a,d)
end
get_targetargs(::Type{<:SST.RegularGridTree}) = (-180.0:30.0:180.0,90.0:-30.0:-90.0)
get_targetkwargs(::Type{<:SST.RegularGridTree}) = (;chunksize=5)
get_targetargs(::Type{<:SST.ISEACircleTree}) = (2,1)
get_targetkwargs(::Type{<:SST.ISEACircleTree}) = (;)
get_targetargs(::Type{<:SST.WebMercatorTree}) = (4,)
get_targetkwargs(::Type{<:SST.WebMercatorTree}) = (;chunksize=4)
get_sourceargs(::Type{<:SST.ISEACircleTree}) = ()
get_sourceargs(::Type{<:SST.RegularGridTree}) = ((:lon,:lat),)


@testset "Reprojection" begin

testindices_3d = [
    (1,1,1),
    (1:2,1,1),
    (1:3,2:4,2),
    (2:4,3:4,3:3),
    (1,3,3:4),
    (:,:,:),
    (2:3,1:2,3:5),
]

testindices_2d = [
    (1,1),
    (2,3),
    (1:3,2:4),
    (1:2,4),
    (3,3:4),
    (:,:),
]

sourcetypes = [SST.RegularGridTree, SST.ISEACircleTree]
targettypes = [SST.RegularGridTree, SST.ISEACircleTree, SST.WebMercatorTree]

for sourcetype in sourcetypes
    for targettype in targettypes
    @testset "$sourcetype to $targettype" begin
        targetargs = get_targetargs(targettype)
        targetkwargs = get_targetkwargs(targettype)
        sourceargs = get_sourceargs(sourcetype)

        a = make_testarray(sourcetype)
        source = SST.ProjectionSource(sourcetype, a, sourceargs... )
        target = SST.ProjectionTarget(targettype, targetargs...; targetkwargs...)

        projar = @test_nowarn SST.LazyProjectedDiskArray(source, target)

        # We first compute the "true" results by reprojecting every point from target to source
        targetcoords_unitsphere = SST.index_to_unitsphere.(LinearIndices(projar), (target.tree,))
        targetcoords_sourcecrs = inv(SST.get_projection(source.tree)).(targetcoords_unitsphere)
        closest_inds = map(targetcoords_sourcecrs) do tc
            map(source.lookups, tc) do lk,t
                DD.selectindices(lk, DD.Near(t))
            end |> CartesianIndex
        end
        closest_vals = a.data[closest_inds]

        testindices = if ndims(projar) == 2
            testindices_2d
        elseif ndims(projar) == 3
            testindices_3d
        else
            error("Only two or three-dimensional arrays are supported. You have given an array with $(ndims(projar)) dimensions.")
        end

        for inds in testindices
            @testset let inds = inds
                @test projar[inds...] == closest_vals[inds...]
            end
        end
        # Now let's try to reproject into a file
        p = tempname()
        output = SST.create_dataset(target, p, arrayname=:data, backend=:zarr)
        SST.reproject!(output, source, target)
        zarrgroup = Zarr.zopen(p)
        @test Array(zarrgroup["data"]) == closest_vals
    end
    end
end
end