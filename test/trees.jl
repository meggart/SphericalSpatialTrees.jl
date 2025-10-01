using SphericalSpatialTrees: SphericalSpatialTrees as SST
using DiskArrays: mockchunks
using DimensionalData

@testset "Build and Use Trees" begin
    lon_range = X(180:-1:-180)
    lat_range = Y(90:-1:-90)
    geo_array = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
    geo_array = mockchunks(geo_array, (128, 128))

    source = SST.ProjectionSource(SST.RegularGridTree, geo_array)
    dggs_resolution = 5
    chunk_length = 2
    target = SST.ProjectionTarget(SST.ISEACircleTree, dggs_resolution, chunk_length)
    dggs_a = SST.LazyProjectedDiskArray(source, target)
    dggs_length = 2^dggs_resolution

    @test size(dggs_a) == (dggs_length, dggs_length, 10)
    @test size(dggs_a[1:2, 1:2, 1]) == (2, 2)
    @test size(dggs_a[1:dggs_length, 1:dggs_length, 1]) == (dggs_length, dggs_length)
end