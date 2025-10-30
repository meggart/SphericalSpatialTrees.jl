import Proj
import CoordinateTransformations: Transformation
struct TransformationChannel{T1,T2} <: Function
    isinv::Bool
    transforms::Channel{Tuple{T1,T2}}
end
Base.inv(t::TransformationChannel) = TransformationChannel(!t.isinv,t.transforms)
allow_threaded_transformation(::TransformationChannel) = false
allow_threaded_transformation(_) = true
function TransformationChannel(f)
    testitem = f()
    c = Channel{typeof(testitem)}(Threads.nthreads())
    for _ in 1:Threads.nthreads()
        put!(c,f())
    end
    return TransformationChannel(false, c)
end
(t::TransformationChannel)(x...) = with_transform(t) do tt
    tt(x...)
end
function with_transform(f,t::TransformationChannel)
    tt = take!(t.transforms)
    _tt = t.isinv ? last(tt) : first(tt)
    try 
        f(_tt)
    finally
        put!(t.transforms,tt)
    end
end
#Generic fallback
with_transform(f,t) = f(t)

function init_threaded_proj_epsg(epsg_code)
    TransformationChannel() do
        (
            Proj.Transformation("OGC:84","EPSG:$epsg_code") ∘ GeographicFromUnitSphere(),
            UnitSphereFromGeographic() ∘ inv(Proj.Transformation("OGC:84","EPSG:$epsg_code"))
        )
    end
end

struct MultiZoneProjection{P} <: Transformation
    projections::Vector{P}
end
(p::MultiZoneProjection)((x,y,zone)) = p.projections[zone]((x,y))

function init_threaded_proj_collection_epsg(epsg_codes)
    TransformationChannel() do
        p = MultiZoneProjection([Proj.Transformation("OGC:84","EPSG:$epsg_code") ∘ GeographicFromUnitSphere() for epsg_code in epsg_codes])
        pinv = MultiZoneProjection([UnitSphereFromGeographic() ∘ Proj.Transformation("EPSG:$epsg_code","OGC:84") for epsg_code in epsg_codes])
        p,pinv
    end
end