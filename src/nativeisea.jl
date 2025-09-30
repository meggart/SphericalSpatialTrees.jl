module NativeISEA
using StaticArrays: SVector, @SVector, @SMatrix
using Rotations: RotX, RotZ
using GeometryOps.UnitSpherical: spherical_distance, UnitSphereFromGeographic, GeographicFromUnitSphere, UnitSphericalPoint, slerp
const USP = UnitSphericalPoint
using LinearAlgebra: cross, dot,norm, normalize
import CoordinateTransformations: Transformation, ∘

export ISEA, ISEA20, InvISEA20, ISEA10, InvISEA10


const NORMAL_OFFSET = -0.7946544722917661229555309283275940420265905883092648017549557750084386449717329
const ITRIANGLE_SIDE_SQ = 0.90450849718747371205114670859140952943007729495144071553386215567
const ITRIANGLE_HEIGHT_SIDE = 1.044436448670983612734708092275980120462245664528624497238884055281095681373414
const MAX_TRIANGLE_HEIGHT = Float64(sind(BigFloat(60)))
const TRIANGLE_AREA = 0.6283185307179586
const TRIANGLE_PRODUCT = 0.7608452130361227
const DISTANCE_THRESHOLD = 0.6408518201709885



function _compute_trans_matrix()
    sqrt3 = sqrt(BigFloat(3))
    tm = @SMatrix([1 -sqrt3/3; sqrt3/2 0.5])
    Float64.(tm), Float64.(inv(tm))
end
const trans_matrix, itrans_matrix = _compute_trans_matrix()


function compute_vertices()
    vertices = Vector{SVector{3,BigFloat}}(undef, 12)
    latring = atand(BigFloat(0.5))

    coslat = cosd(latring)
    sinlat = sind(latring)

    vertices[1] = @SVector [BigFloat(0.0), 0, 1]
    vertices[12] = @SVector [BigFloat(0.0), 0, -1]

    for i = 0:4
        northlon = BigFloat(- 72 * i)
        southlon = northlon - 36
        vertices[2+i] = @SVector [sind(northlon) * coslat, cosd(northlon) * coslat, sinlat]
        vertices[7+i] = @SVector [sind(southlon) * coslat, cosd(southlon) * coslat, -sinlat]
    end

    # orientation  = (BigFloat(31.7174744114611), BigFloat(78.8))
    # r1 = RotX(deg2rad(orientation[1]))
    # r2 = RotZ(deg2rad(orientation[2]))
    # rot = r2 * r1

    #rvertices = (rot,).*vertices
    rvertices = vertices
    return UnitSphericalPoint.(rvertices)
end

function A(a, b, c)
    x = triple_product(a, b, c)
    y = (1 + dot(a, b) + dot(b, c) + dot(c, a))
    if abs(y) < 1e-10
        return 2 * atan(x / y), false
    end
    # @show (1 + dot(a, b) + dot(b, c) + dot(c, a))
    2 * atan(x / y),true 
end
triple_product(a, b, c) = dot(a, cross(b, c))



struct ISEATriangle{T}
    A::UnitSphericalPoint{T}
    B::UnitSphericalPoint{T}
    C::UnitSphericalPoint{T}
    kind::Symbol
    #And store midpoints for faster computations
    M::UnitSphericalPoint{T}
    AB::UnitSphericalPoint{T}
    BC::UnitSphericalPoint{T}
    CA::UnitSphericalPoint{T}
    triangles::Vector{@NamedTuple{A::UnitSphericalPoint{T},B::UnitSphericalPoint{T},C::UnitSphericalPoint{T}}}
end
function ISEATriangle(A, B, C, kind, T)
    AB = USP{T}(slerp(A, B, 0.5))
    BC = USP{T}(slerp(B, C, 0.5))
    CA = USP{T}(slerp(C, A, 0.5))
    M = USP{T}((A + B + C) |> normalize)
    A = USP{T}(A)
    B = USP{T}(B)
    C = USP{T}(C)
    triangles = [
        (A=A, B=M, C=CA),(A=A,B=M,C=CA),(A=C,B=M,C=BC),(A=C,B=CA,C=M),
        (A=B,B=M,C=AB),(A=A,B=AB,C=M),(A=B,B=BC,C=M), (A=B,B=BC,C=M)
    ]
    return ISEATriangle(A,B,C,kind,M,AB,BC,CA,triangles)
end
function Base.show(io::IO, ::MIME"text/plain", tri::ISEATriangle)
    corners = GeographicFromUnitSphere().((tri.A, tri.B, tri.C))
    println(io, "$(tri.kind) triangle with corner coordinates $corners")
end
_base(t::ISEATriangle) = t.C-t.B
_height(t::ISEATriangle) = t.A-(t.B+t.C)/2


#Make a type for the neighboring triangles
struct ISEANeighbor
    i::Int #Index of the neighboring triangle
    type::Int # Type of the connection, 0:left/right, same direction, 1:left/right, opposite direction, 2:up/down
end
function makeneighbors()
    #Make neighbors for north triangles
    northneighbors = map(1:5) do i
        ISEANeighbor(mod1(i-1,5),0),ISEANeighbor(mod1(i+1,5),0),ISEANeighbor(i+15,2)
    end
    southneighbors = map(6:10) do i
        ISEANeighbor(mod1(i+1,5)+5,0),ISEANeighbor(mod1(i-1,5)+5,0),ISEANeighbor(21-i,2)
    end
    northtipneighbors = map(11:15) do i
        ISEANeighbor(mod1(i-1,5)+15,1),ISEANeighbor(i+5,1),ISEANeighbor(21-i,2)
    end
    southtipneighbors = map(16:20) do i
        ISEANeighbor(i-5,1),ISEANeighbor(mod1(i+1,5)+10,1),ISEANeighbor(i-15,2)
    end
    return vcat(northneighbors,southneighbors,northtipneighbors,southtipneighbors)
end

"""
    ISEA{T}

A structure to hold the ISEA triangles and their neighbors.
"""
struct ISEA{T}
    triangles::Vector{ISEATriangle{T}}
    neighbors::Vector{NTuple{3,ISEANeighbor}}
end
function ISEA(T=Float64)
    vertices = compute_vertices()
    northindices = [(1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,6,2)]
    southindices = [(12,11,10),(12,10,9),(12,9,8),(12,8,7),(12,7,11)]
    northtriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :north, T) for (i, j, k) in northindices]
    southtriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :south, T) for (i, j, k) in southindices]
    # Generate the triangles for the equator, we make sure that the triangles
    equatorindices_northtip = [(2,11,7),(3,7,8),(4,8,9),(5,9,10),(6,10,11)]
    equatorindices_southtip = [(7,3,2),(8,4,3),(9,5,4),(10,6,5),(11,2,6)]
    equatorindices = vcat(equatorindices_northtip,equatorindices_southtip)
    equatortriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :middle, T) for (i, j, k) in equatorindices]

    triangles = vcat(northtriangles,southtriangles,equatortriangles)
    neighbors = makeneighbors()
    ISEA(triangles,neighbors)
end
Base.show(io::IO,isea::ISEA) = println(io,"ISEA projection with orientation $(GeographicFromUnitSphere()(isea.triangles[1].A))")


const CISEA = ISEA()

"""
    ISEA20(ISEA=ISEA())

A transformation that transforms coordinates from `UnitSphericalPoint`s to x/y coordinates in one of the 20
triangles defining the ISEA projection. Will by defualt use an ISEA instantiation with orienatation to the north 
pole if not passed another one during construction. 
"""
struct ISEA20{T} <: Transformation
    isea::ISEA{T}
end
struct InvISEA20{T} <: Transformation
    isea::ISEA{T}
end
Base.inv(isea::ISEA20) = InvISEA20(isea.isea)
Base.inv(isea::InvISEA20) = ISEA20(isea.isea)
ISEA20(args...;kwargs...) = ISEA20(ISEA(args...;kwargs...))
InvISEA20(args...;kwargs...) = InvISEA20(ISEA(args...;kwargs...))

function to_single_plane((x, y, itri))
    if 1 <= itri <= 5 #North triangles
        return x + itri - 1, y
    elseif 6 <= itri <= 10 # South triangles
        return 1.5 + mod(4 - itri, 5) - x, -y - MAX_TRIANGLE_HEIGHT
    elseif 11 <= itri <= 15 # North tip equator
        return x + mod(itri + 3, 5) + 0.5, y - MAX_TRIANGLE_HEIGHT
    else
        return 1 + (itri - 16) - x, -y
    end
end


function fast_triangle_distance(t::ISEATriangle, p)
    if norm(t.M-p) > DISTANCE_THRESHOLD
        return Inf
    else
        return triangle_distance(t.A,t.B,t.C,p)
    end
end

(isea::ISEA20)(latlon::NTuple{2}) = isea(UnitSphericalPoint(latlon))
function _transform_isea(isea::ISEA,p::UnitSphericalPoint)
    grid = isea
    t = isea.triangles[1]
    d = fast_triangle_distance(t,p)
    iszero(d) && return (first(transform_point(p,t))..., 1)
    current_min = d, 1
    for i in 2:20
        t = grid.triangles[i]
        d = fast_triangle_distance(t,p)
        iszero(d) && return (first(transform_point(p,t))..., i)
        if d < first(current_min)
            current_min = d, i
        end
    end
    #Nothing was found, so lets use the triangle with the smallest distance
    return (first(transform_point(p,grid.triangles[last(current_min)]))..., last(current_min))
end
(isea::ISEA20)(p::UnitSphericalPoint) = _transform_isea(isea.isea,p)

@inline function transform_bary(v0,v1,v2,v)
    #trip = triple_product(v0,v1,v2)
    #p1 = trip*v - triple_product(v,v1,v2)*v0
    p1 = cross(cross(v0,v),cross(v1,v2))
    all(iszero, p1) || (p1 = normalize(p1))
    pstable =  dot(p1,v0) < 0
    h = sqrt((1 - dot(v0, v)) / (1 - dot(v0, p1)))
    a_part, stable = A(v0,v1,p1)
    a_whole, wstable = A(v0,v1,v2)
    β2 = h * a_part / a_whole
    β0 = 1 - h
    β1 = h-β2
    @SVector([β0, β1, β2]), (pstable && wstable && stable)
end
function bary_to_xy(β,rect)
    [rect.A rect.B rect.C] * β
end

"Coordinates of the triangles in the 2d plane"
function _compute_planecoords()
    M = @SVector([0.5, sqrt(3) / 6])
    A = @SVector([0.5, sqrt(3) / 2])
    B = @SVector([0.0, 0.0])
    C = @SVector([1.0, 0.0])
    AB = @SVector([0.25, sqrt(3) / 4])
    BC = @SVector([0.5, 0.0])
    CA = @SVector([0.75, sqrt(3) / 4])
    triangles = (
        (A=A,B=M,C=CA),(A=A,B=M,C=CA),(A=C,B=M,C=BC),(A=C,B=CA,C=M),
        (A=B,B=M,C=AB),(A=A,B=AB,C=M),(A=B,B=BC,C=M),(A=B,B=BC,C=M)
    )
    (;M,A,B,C,AB,BC,CA,triangles)
end
const PlaneCoordinates = _compute_planecoords()

"""
   Determines the rectangular sub-triangle from a given UnitSphericalPoint
"""
function find_subtriangle(t, v)
    b_c = dist_rhs(t.A, t.BC, v) < 0
    c_a = dist_rhs(t.B, t.CA, v) < 0
    a_b = dist_rhs(t.C, t.AB, v) < 0
    ((b_c << 2) | (c_a << 1) | a_b) + 1
end

function transform_point(v, t)
    itri = find_subtriangle(t,v)
    trirect = t.triangles[itri]
    planerect = PlaneCoordinates.triangles[itri]   
    β, stable = transform_bary(trirect..., v)
    bary_to_xy(β,planerect),stable, itri
end

function intriangle(p,tri::ISEATriangle) 
    iszero(triangle_distance(tri.A,tri.B,tri.C,p))
end

"Projected signed distance of a point v is to the right of the line between a and b"
dist_rhs(a::UnitSphericalPoint,b::UnitSphericalPoint,v::UnitSphericalPoint) = dot(cross(a,b),v)

function dist_rhs(a::SVector{2},b::SVector{2},v::SVector{2})
    c = b-a
    c[1]*(v[2]-a[2]) - c[2]*(v[1]-a[1])
end


function triangle_distance(a,b,c,v)
    d_ab = dist_rhs(b,a,v)
    d_bc = dist_rhs(c,b,v)
    d_ca = dist_rhs(a,c,v)
    max(d_ab,d_bc,d_ca,0.0)
end


# function itransform_point(x1,x2,tri::ISEATriangle)
#     base=_base(tri)
#     height=_height(tri)
#     pointonsurface = tri.B + x1 * base + x2 * height / MAX_TRIANGLE_HEIGHT
#     p = pointonsurface/norm(pointonsurface)
#     return GeographicFromUnitSphere()(p)
# end
function itransform_point(t, x, y)
    itri = find_subtriangle(PlaneCoordinates,@SVector([x,y]))
    planetri = PlaneCoordinates.triangles[itri]
    trimat = vcat(@SMatrix([1.0 1.0 1.0]),[planetri.A planetri.B planetri.C])
    β = trimat \ @SVector([1.0,x,y])
    β0, _, β2 = β
    h = 1 - β0
    t = t.triangles[itri]
    iszero(h) && return t.A
    q = if β2 != 0.0
        a = β2 / h * first(A(t.A,t.B,t.C))
        S = sin(a)
        C = 1 - cos(a)
        #@show a,S,C
        f = S * triple_product(t.A,t.B,t.C) + C * (dot(t.A, t.B) * dot(t.B, t.C) - dot(t.C, t.A))
        g = C * sqrt(1 - dot(t.B, t.C)^2) * (1 + dot(t.A, t.B))
        #@show f,g
        2 / acos(dot(t.B, t.C)) * atan(g / f)
    else
        0.0
    end
    p = slerp(t.B, t.C, q)
    T = acos(1 + h^2 * (dot(t.A, p) - 1)) / acos(dot(t.A, p))
    slerp(t.A, p, T)
end
(isea::InvISEA20)((x1, x2, n)) = itransform_point(isea.isea.triangles[n], x1, x2)


const northpairs = [(i,i+15) for i in 1:5]
const southpairs = [(i,21-i) for i in 11:15]
const diamond2tri = vcat(northpairs,southpairs[[2,3,4,5,1]])
const tri2diamond = [(1,1),(2,1),(3,1),(4,1),(5,1),
               (9,2),(8,2),(7,2),(6,2),(10,2),
               (10,1),(6,1),(7,1),(8,1),(9,1),
               (1,2),(2,2),(3,2),(4,2),(5,2)]

struct ISEATrianglesToDiamond <: Transformation end
struct ISEADiamondToTriangles <: Transformation end
Base.inv(::ISEATrianglesToDiamond) = ISEADiamondToTriangles()
Base.inv(::ISEADiamondToTriangles) = ISEATrianglesToDiamond()

function (::ISEATrianglesToDiamond)((x1,x2,n))
    i,j = tri2diamond[n]
    if j == 2
        #We are in the south diamond, so we need to reverse the x2 coordinate
        x1 = 1-x1
        x2 = -x2
    end
    return x1,x2,i
end
function (::ISEADiamondToTriangles)((x1,x2,i))
    i1,i2 = diamond2tri[i]
    i,x1,x2 = x2 < 0 ? (i2,1-x1,-x2) : (i1,x1,x2)
    x1,x2,i
end

ISEA10(isea::ISEA) = ISEATrianglesToDiamond() ∘ ISEA20(isea)
ISEA10(args...;kwargs...) = ISEATrianglesToDiamond() ∘ ISEA20(args...;kwargs...)

struct RotateISEA <: Transformation end
struct InvRotateISEA <: Transformation end
Base.inv(::RotateISEA) = InvRotateISEA()
Base.inv(::InvRotateISEA) = RotateISEA()

function (::RotateISEA)((x,y,i))
    x,y = trans_matrix * @SVector([x,y])
    x,y/sqrt(3)*2,i
end
function (::InvRotateISEA)((x,y,i))
    x1, x2 = NativeISEA.itrans_matrix * @SVector([x, y * sqrt(3) / 2])
    x1,x2,i
end


struct ISEARectToDiamond <: Transformation end
struct ISEADiamondToRect <: Transformation end
Base.inv(::ISEARectToDiamond) = ISEADiamondToRect()
Base.inv(::ISEADiamondToRect) = ISEARectToDiamond()


function (::ISEARectToDiamond)((x1,x2,i) )
    isdown = dist_rhs(@SVector([0.5,-sqrt(3)/2]),@SVector([1.0,0.0]),@SVector([x1,x2])) < 0
    i = mod1(i+1,5)
    if isdown
        x1 = x1-0.5
        x2 = x2+sqrt(3)/2
        i = i + 5
    end
    x1,x2,i
end
_diamond2rect(i) = mod1(i-1,5),(i-1)÷5+1
function (::ISEADiamondToRect)((x,y,i))
    i,lr = _diamond2rect(i)
    if lr == 2
        x = 0.5+x
        y = -sqrt(3)/2 + y
    end
    x,y,i
end

ISEA5(isea::ISEA) = ISEADiamondToRect() ∘ ISEATrianglesToDiamond() ∘ ISEA20(isea)
ISEA5(args...;kwargs...) = ISEADiamondToRect() ∘ ISEATrianglesToDiamond() ∘ ISEA20(args...;kwargs...)


struct PickPlane <: Transformation
    i::Int
end
(p::PickPlane)(coords) = (coords...,p.i)



# struct ISEA5Cell{T} <: Function
#     isea::ISEA{T}
#     resolution::Int
# end
# ISEA5Cell(resolution::Int,args...;kwargs...) = ISEA5Cell(ISEA(args...;kwargs...),resolution)
# function (isea::ISEA5Cell)(p::UnitSphericalPoint)
#     n,x1,x2 = ISEA5(isea.isea)(p)
#     x1rot,x2rot = trans_matrix * @SVector([x1,x2])
#     i_cell = floor(Int,x1rot * 2^isea.resolution)
#     j_cell = floor(Int,x2rot / sqrt(3)*2 * 2^isea.resolution)
#     i_cell = clamp(i_cell, 0, 2 * 2^isea.resolution - 1)
#     j_cell = clamp(j_cell, 0, 2^isea.resolution - 1)
#     Cell(n,i_cell,j_cell,isea.resolution)
# end
# (isea::ISEA5Cell)(latlon::NTuple{2}) = isea(UnitSphericalPoint(latlon))

# struct InvISEA5Cell{T} <: Function
#     isea::ISEA{T}
#     resolution::Int
# end
# Base.inv(isea::ISEA5Cell) = InvISEA5Cell(isea.isea,isea.resolution)
# Base.inv(isea::InvISEA5Cell) = ISEA5Cell(isea.isea,isea.resolution)
# function (isea::InvISEA5Cell)(c::Cell) 
#     c.resolution == isea.resolution || throw(ArgumentError("Trying to convert with wrong resolution"))
#     x = c.i / 2^isea.resolution
#     y = c.j / 2^isea.resolution * sqrt(3)/2
#     x,y = itrans_matrix * @SVector([x,y])
#     return InvISEA5(isea.isea)((c.n,x,y))
# end



end