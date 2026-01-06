module SphericalSpatialTreesGeoMakieExt

import SphericalSpatialTrees as SST
import GeoMakie

GeoMakie.Makie.args_preferred_axis(::Type{<: GeoMakie.Makie.Poly}, ::SST.LazyProjectedDiskArray) = GeoMakie.GlobeAxis

end