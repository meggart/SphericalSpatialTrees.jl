import CoordinateTransformations: Transformation, ∘
const R2D = 180 / pi
const RE = 6378137.0
#Webmercator projection code copied from MapTiles.jl
#https://github.com/JuliaGeo/MapTiles.jl/blob/1a838193228b45579130bba940d68d03397b97a5/src/tiles.jl#L37-L63
struct WebMercatorToLonLat <: Transformation end
function (::WebMercatorToLonLat)((x,y))
    lng = x * R2D / RE
    lat = ((pi / 2) - 2.0 * atan(exp(-y / RE))) * R2D
    return lng, lat
end
struct LonLatToWebMercator <: Transformation end
function (::LonLatToWebMercator)((lng,lat))
    x = RE * deg2rad(lng)
    y = if lat <= -90
        -Inf
    elseif lat >= 90
        Inf
    else
        RE * log(tan((pi / 4) + (0.5 * deg2rad(lat))))
    end
    return x, y
end
Base.inv(::LonLatToWebMercator) = WebMercatorToLonLat()
Base.inv(::WebMercatorToLonLat) = LonLatToWebMercator()


const halfsquarelength = 2.0037508342789244e7
struct WebMercatorTree end
function WebMercatorTree(max_level)
    r = range(-halfsquarelength,halfsquarelength,length=2^max_level+1)
    p = UnitSphereFromGeographic() ∘ WebMercatorToLonLat()
    RegularGridTree(r,r,p)
end

function ProjectionTarget(::Type{<:WebMercatorTree},target_resolution)
    tree = WebMercatorTree(target_resolution)
    chunktree = WebMercatorTree(target_resolution-8)
    ProjectionTarget(tree,chunktree)
end
Base.ndims(::Type{WebMercatorTree}) = 2