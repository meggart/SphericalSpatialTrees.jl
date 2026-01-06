import SphericalSpatialTrees as SST
import GeometryOps as GO, GeoInterface as GI
import DiskArrays

using GeometryOps.UnitSpherical: GeographicFromUnitSphere, spherical_distance, 
                            UnitSphereFromGeographic, SphericalCap

using Rasters, RasterDataSources, ArchGDAL, Zarr # data sources
using GeoMakie, GLMakie # visualization

get!(ENV, "RASTERDATASOURCES_PATH", mkpath(joinpath(first(Base.DEPOT_PATH), "artifacts", "RasterDataSources"))) # hide
# nothing # hide

# Now that we have loaded all of these packages, let's start talking about functionality.
#
# First, let's use the top-level interface to create a lazily regridded array, from a
# lat-long source to a DGGS target.
#
# This dataset is a bioclimatic dataset, which is sufficiently small
# that you can do this quickly on your machine.

ras = Raster(EarthEnv{HabitatHeterogeneity}, :cv) |> r -> replace_missing(r, NaN) .|> log

# Let's give this dataset some chunks.  It's 1728x696 so 
# 75x75 chunks give it ~ 22x11 chunks.

ras_chunked = DiskArrays.mockchunks(ras, (75, 75))

# Now, we need to define the source and target trees.  The dataset itself is
# on a plane (equirectangular projection), so we can use the `RegularGridTree`
# which is an implicit quadtree.

source = SST.ProjectionSource(SST.RegularGridTree, ras_chunked);

# We want to go to a DGGS target, so we can use the `ISEACircleTree` 
# (based on the Snyder equal area projection).  This is the only tree so far
# but more can be added as we go forward.

target = SST.ProjectionTarget(SST.ISEACircleTree, 8, 2)

# As a sidenote, this is what actually happens internally.  Using the spatial
# tree interface, we compute the map of chunks in the source dataset that each
# chunk in the target dataset requires to be materialized.

SST.compute_connected_chunks(source, target)

# Now, we can create the lazy projected disk array.

a = SST.LazyProjectedDiskArray(source, target)

# You can now access some data, using this as you would use any array:

a[1:10,1:10,1]

# You can also materialize the array fully into memory (if the array is small enough):

ac = collect(a);

# Now let's visualize this!  GLMakie has no problem visualizing this many polygons:
# First, let's get every polygon:

polys = SST.index_to_polygon_lonlat.(vec(collect(eachindex(a))), (target.tree,));

# Then we can just assign the correct color to the correct polygon,
# and plot on GeoMakie's `GlobeAxis`.

fig, ax, p1 = poly(a; axis = (; type = GlobeAxis))
p2 = lines!(ax, a; color = :black, transparency = true, linewidth = 0.075)
p3 = meshimage!(ax, -180..180, -90..90, fill(colorant"white", 2, 2); zlevel = -100_000)
fig
    
display(GLMakie.Screen(), meshimage(-180..180, -90..90, reorder(ras, Y => Rasters.ForwardOrdered()); axis = (; type = GlobeAxis)))

# To give a better idea of what this looks like, let's use a lower level (lower resolution) DGGS:

target = SST.ProjectionTarget(SST.ISEACircleTree, 3, 2)
a = SST.LazyProjectedDiskArray(source, target)
ac = mapreduce((i,j)->cat(i,j,dims=3),1:10) do n
    a[:,:,n]
end  
fig, ax, plt = poly(a; strokewidth = 1, strokecolor = :black, axis = (; type = GlobeAxis));
meshimage!(ax, -180..180, -90..90, fill(colorant"white", 2, 2); zlevel = -100_000) # background plot
fig

fig2, ax2, plt2 = meshimage(-180..180, -90..90, reorder(ras, Y => Rasters.ForwardOrdered()); zlevel = 100_000, alpha = 0.5, transparency = true, axis = (; type = GlobeAxis, show_axis = false))
poly!(ax2, a; strokewidth = 1, strokecolor = :black)
meshimage!(ax2, -180..180, -90..90, fill(colorant"white", 2, 2); zlevel = -100_000) # background plot

# And that's the basics!
# Now we can go a bit into the weeds of capabilities and how this works.

mycoord = (11.0, 53.0)
rad = 1e-6

iseat = SST.ISEACircleTree(22)
mycircle = SphericalCap(UnitSphereFromGeographic()(mycoord), rad)

# Query

r = GO.query(iseat, mycircle)
points = SST.index_to_lonlat.(r, (iseat,))
f, a, p = scatter(points; color = Makie.Cycled(2), axis = (; type = GlobeAxis, show_axis = false));
lines!(a, GeoMakie.coastlines())
f

#
poly1 = SST.index_to_polygon_lonlat(1, iseat)

poly!(a, poly1; strokecolor = :red, strokewidth = 1)
f

#
lo,up = CartesianIndices(SST.gridsize(iseat))[r] |> extrema
bb_isea = map(Colon(), lo.I, up.I)


ras = Raster(WorldClim{BioClim}, 5) |> r -> replace_missing(r, NaN)
c = DiskArrays.mockchunks(ras, DiskArrays.GridChunks(size(ras), (100, 100)))
lccst = SST.RegularGridTree(c)
r2 = GO.query(SST.rootnode(lccst), mycircle)
SST.index_to_lonlat.(r2, (lccst,))

lo, up = CartesianIndices(size(c)[1:2])[r2] |> extrema
bb_lccs = map(Colon(), lo.I, up.I)



d, i = SST.find_nearest(iseat, (11.0, 53.0))
SST.index_to_lonlat(i, iseat)

d, i = SST.find_nearest(lccst, (11.0, 53.0))
SST.index_to_lonlat(i, lccst)



source = SST.ProjectionSource(SST.RegularGridTree, c);
target = SST.ProjectionTarget(SST.ISEACircleTree, 3, 2);


SST.compute_connected_chunks(source, target)


a = SST.LazyProjectedDiskArray(source, target)

# How do I get a polygon from a cartesian index?

SST.index_to_polygon_lonlat(CartesianIndex(1,1,1), target.tree)

# Let's get every polygon...

polys = SST.index_to_polygon_lonlat.(eachindex(a), (target.tree,));

# This is 

length(polys)

#
a[:, 1, :]

# Now we can materialize the lazy array...

@time ac = collect(a);

# and GLMakie has no problem visualizing this!

poly(vec(polys); color = vec(ac), strokewidth = 1, strokecolor = :black, axis = (; type = GlobeAxis, show_axis = false))
meshimage!(-180..180, -90..90, fill(colorant"white", 2, 2); zlevel = -100_000)
current_figure() 

# heatmap(c.lccs_class.data[bb_lccs...][:, :, 1])



# ds = SST.create_dataset(target, "./output.zarr/", arrayname=:lccs_class)

# SST.reproject!(ds.lccs_class,source,target)
