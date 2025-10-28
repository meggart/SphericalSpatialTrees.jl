import Proj
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
function init_threaded_proj_epsg(epsg_code)
    TransformationChannel() do
        (
            Proj.Transformation("OGC:84","EPSG:$epsg_code") ∘ GeographicFromUnitSphere(),
            UnitSphereFromGeographic() ∘ inv(Proj.Transformation("OGC:84","EPSG:$epsg_code"))
        )
    end
end