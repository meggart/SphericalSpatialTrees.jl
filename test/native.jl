using GeometryOps.UnitSpherical: UnitSphereFromGeographic, spherical_distance, GeographicFromUnitSphere

using SphericalSpatialTrees.NativeISEA: ISEA5, ISEA10, ISEA20
@testset "Native implementation" begin
    for isea in [ISEA5(), ISEA10(), ISEA20()]
        lon_range = -180.0:1:180.0
        lat_range = 90.0:-1:-90.0
        lonlatcoords = tuple.(lon_range, lat_range')
        res = isea.(lonlatcoords)
        back = inv(isea).(res)
        dist = spherical_distance.(UnitSphereFromGeographic().(lonlatcoords), back)
        @test maximum(dist) < 1e-7
    end
end
