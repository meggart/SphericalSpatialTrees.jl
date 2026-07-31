import Proj
import CoordinateTransformations: Transformation
# Mutable wrapper around a channel that has a finalizer for projection contexts
_tuple_types(::Type{Tuple{T1,T2}}) where {T1,T2} = T1, T2
mutable struct _TransformChannel{T1,T2}
    transforms::Channel{Tuple{T1,T2}}
    contexts::Vector{Any}
    function _TransformChannel(transforms::Channel, contexts)
        ct = eltype(transforms)
        T1, T2 = _tuple_types(ct)
        x = new{T1,T2}(transforms, contexts)
        finalizer(x) do t
            foreach(Proj.proj_context_destroy, t.contexts)
        end
        x
    end
end

struct TransformationChannel{T1,T2} <: Function
    isinv::Bool
    transforms::_TransformChannel{T1,T2}
end
Base.inv(t::TransformationChannel) = TransformationChannel(!t.isinv, t.transforms)
allow_threaded_transformation(::TransformationChannel) = true
allow_threaded_transformation(_) = true
function TransformationChannel(f)
    testitem = f(C_NULL)
    c = Channel{typeof(testitem)}(Threads.nthreads())
    contexts = Any[]
    for _ in 1:Threads.nthreads()
        # We create one context per thread
        ctx = Proj.proj_context_create()
        put!(c, f(ctx))
        push!(contexts, ctx)
    end
    TransformationChannel(false, _TransformChannel(c, contexts))
end
(t::TransformationChannel)(x...) = with_transform(t) do tt
    tt(x...)
end
function with_transform(f,t::TransformationChannel)
    tt = take!(t.transforms.transforms)
    _tt = t.isinv ? last(tt) : first(tt)
    try 
        f(_tt)
    finally
        put!(t.transforms,tt)
    end
end
#Generic fallback
with_transform(f,t) = f(t)

function init_threaded_proj(crs)
    TransformationChannel() do ctx
        (
            Proj.Transformation("OGC:84", crs, always_xy=true, ctx=ctx) ∘ GeographicFromUnitSphere(),
            UnitSphereFromGeographic() ∘ inv(Proj.Transformation("OGC:84", crs, ctx=ctx))
        )
    end
end

struct MultiZoneProjection{P} <: Transformation
    projections::Vector{P}
end
(p::MultiZoneProjection)((x,y,zone)) = p.projections[zone]((x,y))

function init_threaded_proj_collection(crss)
    TransformationChannel() do ctx
        p = MultiZoneProjection([Proj.Transformation("OGC:84", crs, ctx=ctx) ∘ GeographicFromUnitSphere() for crs in crss])
        pinv = MultiZoneProjection([UnitSphereFromGeographic() ∘ Proj.Transformation(crs, "OGC:84", ctx=ctx) for crs in crss])
        p,pinv
    end
end