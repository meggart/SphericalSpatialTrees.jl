import Proj

const halfsquarelength = 2.0037508342789244e7
struct WebMercatorTree end
function WebMercatorTree(max_level)
    r = range(-halfsquarelength,halfsquarelength,length=2^max_level+1)
    p = UnitSphereFromGeographic() ∘ Proj.Transformation("+proj=webmerc +datum=WGS84","OGC:84")
    RegularGridTree(r,r,p)
end

function ProjectionTarget(::Type{<:WebMercatorTree},target_resolution)
    tree = WebMercatorTree(target_resolution)
    chunktree = WebMercatorTree(target_resolution-8)
    ProjectionTarget(tree,chunktree)
end
Base.ndims(::Type{WebMercatorTree}) = 2