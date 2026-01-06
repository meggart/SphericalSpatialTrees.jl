module SphericalSpatialTreesMakieExt
import SphericalSpatialTrees as SST
import Makie

Makie.convert_arguments(::Type{<:Makie.Poly}, a::SST.LazyProjectedDiskArray) = (a,)

function Makie.plot!(p::Makie.Poly{<: Tuple{<: SST.LazyProjectedDiskArray}})
    Makie.map!(p, [:polygon], [:polygons, :polygons_colors]) do a
        idxs = vec(collect(eachindex(a)))
        polys = SST.index_to_polygon_lonlat.(idxs, (a.target.tree,))
        return (polys, a[idxs])
    end
    Makie.poly!(p, p.attributes, p[:polygons]; color = p[:polygons_colors])
end

function Makie.convert_arguments(::Makie.PointBased, a::SST.LazyProjectedDiskArray)
    Makie.convert_arguments(Makie.PointBased(), SST.index_to_polygon_lonlat.(eachindex(a), (a.target.tree,)))
end

end