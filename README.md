# SphericalSpatialTrees.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://meggart.github.io/SphericalSpatialTrees.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://meggart.github.io/SphericalSpatialTrees.jl/dev/)
[![Build Status](https://github.com/meggart/SphericalSpatialTrees.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/meggart/SphericalSpatialTrees.jl/actions/workflows/CI.yml?query=branch%3Amain)

SphericalSpatialTrees.jl is a Julia package for regridding and querying geospatial data between different projections, e.g., between rasters in a geographical projection and a Discrete Global Grid System.  It provides spatial tree data structures and projection utilities for this purpose.

Required chunks from the other array are identified using quad tree search lazily on request.
Chunk intersection is tested based on the spherical caps of their chunk boundaries.
The regridding itself is currently done using nearest neighbor search, though other methods (area-conserving, etc.) could also be used. 

## Get Started

Install the package:

```julia
using Pkg
Pkg.add(url="https://github.com/meggart/SphericalSpatialTrees.jl")
```

Create example raster data in equirectangular projection:

```julia
import DiskArrays
using DimensionalData

lons = X(180:-1:-180)
lats = Y(90:-1:-90)
unchunked_geo_array = [exp(cosd(lon)) + 3 * (lat / 90) for lon in lons, lat in lats] # create a 'raster' dataset
geo_array = DiskArrays.mockchunks(unchunked_geo_array, (128, 128)) # create 'fake' chunks
```

Lazy regridding to a Discrete Global Grid System:

```julia
import SphericalSpatialTrees as SST
dggs_resolution = 5
chunk_length = 2
source = SST.ProjectionSource(SST.RegularGridTree, geo_array)
target = SST.ProjectionTarget(SST.ISEACircleTree, dggs_resolution, chunk_length)
dggs_array = SST.LazyProjectedDiskArray(source, target)

# 32×32×10 LazyProjectedDiskArray{Float64}
```

Start reprojection:

```julia
dggs_array[1:8,1:8,1]

# 8×8 Matrix{Float64}:
#  1.8827   1.94937  2.01603  2.0827   2.14937  2.21603  2.2827   2.34937
#  1.79904  1.84901  1.94901  1.99929  2.06595  2.13262  2.19929  2.26595
#  1.73262  1.7832   1.84987  1.90075  1.96742  2.03408  2.10075  2.15193
#  1.64987  1.70075  1.75193  1.8186   1.90341  1.97008  2.02185  2.08852
#  1.56742  1.6186   1.67008  1.75519  1.80726  1.87393  1.92629  1.99296
#  1.5186   1.57008  1.60726  1.67393  1.72629  1.77895  1.87895  1.93189
#  1.43674  1.47393  1.52629  1.57895  1.63189  1.73189  1.78512  1.83863
#  1.35519  1.39296  1.44561  1.49856  1.58512  1.63863  1.69242  1.75909
```
