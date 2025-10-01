import SphericalSpatialTrees as SST
import GeometryOps.UnitSpherical: GeographicFromUnitSphere, UnitSphereFromGeographic, SphericalCap, spherical_distance
import GeometryOps.SpatialTreeInterface as STI
using Test

wmt = SST.WebMercatorTree(20)

mycoord = (11.0, 53.0)
rad = 0.000002

mycircle = SphericalCap(UnitSphereFromGeographic()(mycoord), rad)

r = STI.query(SST.rootnode(wmt), mycircle)
@test r == [741342739752, 741342739753, 741343788328, 741343788329, 741344836904, 741344836905]
lonlats = SST.index_to_lonlat.(r,(wmt,))
@test lonlats == [(10.999889373779299, 52.99980605969741),(11.000232696533207, 52.99980605969741),(10.999889373779299, 53.00001267692221),
 (11.000232696533203, 53.00001267692221),(10.999889373779299, 53.00021929315826),(11.000232696533205, 53.00021929315826)]

target = SST.ProjectionTarget(SST.WebMercatorTree,20)

@test target.tree isa SST.RegularGridTree
@test target.tree.x == -2.0037508342789244e7:38.21851414258813:2.0037508342789244e7
@test target.tree.y == -2.0037508342789244e7:38.21851414258813:2.0037508342789244e7
@test target.tree.trans((-2.0037508342789244e7,-2.0037508342789244e7)) |> GeographicFromUnitSphere() == (-180.0, -85.05112877980663)
