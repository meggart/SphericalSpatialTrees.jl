using Test
import SphericalSpatialTrees as SST
import DimensionalData as DD
import DiskArrays

# To exercise the sequential projection path in readblock! we need a fine source
# grid whose chunks are numerous enough that a single target chunk connects to
# >=100 source chunks (threshold in LazyProjection.jl:258).
#
# Source: RegularGridTree 720x360, chunked 10x10  => 72x36 = 2592 source chunks
# Target: RegularGridTree  18x9,  chunked  3x3  =>   6x3 = 18 target chunks
# Each target chunk spans ~40x40 source cells, triggering 200-480 connected chunks.

@testset "LazyProjection — sequential path" begin
    a = DiskArrays.mockchunks(rand(Float32, 720, 360), (10, 10))
    d = DD.Dim{:lon}(-179.75:0.5:179.75), DD.Dim{:lat}(89.75:-0.5:-89.75)
    da = DD.DimArray(a, d)

    source = SST.ProjectionSource(SST.RegularGridTree, da, (:lon, :lat))
    target = SST.ProjectionTarget(
        SST.RegularGridTree,
        -180.0:20.0:180.0,
        90.0:-20.0:-90.0;
        chunksize=3,
    )

    @test SST.nleaf(source.chunktree) == 2592
    @test SST.nleaf(target.chunktree) == 18

    lazy = SST.LazyProjectedDiskArray(source, target)

    # Reference values via direct nearest-neighbor lookup
    targetcoords_unitsphere = SST.index_to_unitsphere.(
        LinearIndices(lazy), (target.tree,)
    )
    targetcoords_sourcecrs = inv(SST.get_projection(source.tree)).(
        targetcoords_unitsphere
    )
    closest_inds = map(targetcoords_sourcecrs) do tc
        map(source.lookups, tc) do lk, t
            DD.selectindices(lk, DD.Near(t))
        end |> CartesianIndex
    end
    expected = a.parent[closest_inds]

    # Sequential path: block reads (each triggers precompute_sequential_weights)
    @test lazy[1:3, 1:3] == expected[1:3, 1:3]
    @test lazy[4:6, 1:3] == expected[4:6, 1:3]
    @test lazy[7:9, 4:6] == expected[7:9, 4:6]
    @test lazy[10:12, 4:6] == expected[10:12, 4:6]
    @test lazy[13:15, 7:9] == expected[13:15, 7:9]
    @test lazy[16:18, 7:9] == expected[16:18, 7:9]
    @test lazy[:, :] == expected
end

@testset "LazyProjection — index buffer growth" begin
    a = DiskArrays.mockchunks(rand(Float32, 720, 360), (10, 10))
    d = DD.Dim{:lon}(-179.75:0.5:179.75), DD.Dim{:lat}(89.75:-0.5:-89.75)
    da = DD.DimArray(a, d)

    source = SST.ProjectionSource(SST.RegularGridTree, da, (:lon, :lat))
    target = SST.ProjectionTarget(
        SST.RegularGridTree,
        -180.0:20.0:180.0,
        90.0:-20.0:-90.0;
        chunksize=3,
    )

    # Default buffer size is 100, but each target chunk connects to 200+ source
    # chunks. The buffer growth path in precompute_sequential_weights lines 39-43
    # is exercised when the readblock! index_arraybuffer is too small.
    buf = SST.make_indexbuffer(source.tree, target.tree)
    @test length(buf) == 100

    chunks = SST.compute_connected_chunks(source, target, (1:3, 1:3))
    @test length(chunks) > 100

    # Full read via DimArray interface (sequential for every chunk)
    lazy = SST.LazyProjectedDiskArray(source, target)
    @test_nowarn SST.DD.DimArray(lazy, SST.DD.dims(target.tree))
end

@testset "LazyProjection — utilities" begin
    a = DiskArrays.mockchunks(rand(Float32, 720, 360), (10, 10))
    d = DD.Dim{:lon}(-179.75:0.5:179.75), DD.Dim{:lat}(89.75:-0.5:-89.75)
    da = DD.DimArray(a, d)

    source = SST.ProjectionSource(SST.RegularGridTree, da, (:lon, :lat))
    target = SST.ProjectionTarget(
        SST.RegularGridTree,
        -180.0:20.0:180.0,
        90.0:-20.0:-90.0;
        chunksize=3,
    )

    # ProjectionTarget show method
    @test repr(MIME("text/plain"), target) isa String

    # 2-arg compute_connected_chunks (full matrix)
    connected = SST.compute_connected_chunks(source, target)
    @test length(connected) == SST.nleaf(target.chunktree)
    @test all(c -> length(c) > 0, connected)

    # DD.DimArray shortcut (source -> target)
    da = SST.DD.DimArray(source, target)
    @test da isa SST.DD.DimArray
    @test size(da) == SST.gridsize(target.tree)
end
