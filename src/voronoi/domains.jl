# Domains that Voronoi diagrams are restricted to.

abstract type AbstractDomain{D} end

"""
    PolygonDomain(vertices)

Simple (possibly non-convex) polygonal domain in 2D. The vertices are stored
counter-clockwise. Boundary faces on edge `k` (from vertex `k` to `k+1`) get
boundary tag `k`.
"""
struct PolygonDomain <: AbstractDomain{2}
    vertices::Vector{Vec2}
    function PolygonDomain(vertices)
        v = to_points(vertices, 2)
        if length(v) > 1 && v[1] ≈ v[end]
            pop!(v)
        end
        length(v) >= 3 || throw(ArgumentError("Polygon must have at least 3 vertices"))
        if signed_area(v) < 0
            reverse!(v)
        end
        return new(v)
    end
end

"""
    ConvexPolyhedronDomain(normals, offsets)

Convex polyhedral domain `{x : dot(normals[k], x) <= offsets[k]}` in 3D.
Boundary faces on facet `k` get boundary tag `k`. See [`box_domain`](@ref)
and [`convex_hull_domain`](@ref).
"""
struct ConvexPolyhedronDomain <: AbstractDomain{3}
    normals::Vector{Vec3}
    offsets::Vector{Float64}
    polytope::ClipPolyhedron
end

function ConvexPolyhedronDomain(normals::AbstractVector, offsets::AbstractVector)
    unit_normals = [Vec3(n) / norm(n) for n in normals]
    unit_offsets = [Float64(o) / norm(n) for (n, o) in zip(normals, offsets)]
    P = polytope_from_planes(unit_normals, unit_offsets)
    return ConvexPolyhedronDomain(unit_normals, unit_offsets, P)
end

"""
    box_domain(lo, hi)

Axis-aligned box `[lo, hi]` in 2D or 3D. In 2D the boundary tags are 1:
`y = lo`, 2: `x = hi`, 3: `y = hi`, 4: `x = lo`. In 3D they are 1/2: `x =
lo/hi`, 3/4: `y = lo/hi`, 5/6: `z = lo/hi`.
"""
function box_domain(lo, hi)
    D = length(lo)
    length(hi) == D || throw(ArgumentError("lo and hi must have the same length"))
    if D == 2
        return PolygonDomain([Vec2(lo[1], lo[2]), Vec2(hi[1], lo[2]), Vec2(hi[1], hi[2]), Vec2(lo[1], hi[2])])
    elseif D == 3
        normals = [Vec3(-1, 0, 0), Vec3(1, 0, 0), Vec3(0, -1, 0), Vec3(0, 1, 0), Vec3(0, 0, -1), Vec3(0, 0, 1)]
        offsets = [-lo[1], hi[1], -lo[2], hi[2], -lo[3], hi[3]]
        return ConvexPolyhedronDomain(normals, offsets)
    else
        throw(ArgumentError("Only 2D and 3D domains are supported"))
    end
end

"""
    convex_hull_domain(points)

Convex domain given as the convex hull of `points` (3D), or the convex hull
polygon in 2D.
"""
function convex_hull_domain(points)
    pts = first(points) isa Real ? throw(ArgumentError("Expected a collection of points")) : points
    D = points isa AbstractMatrix ? size(points, 2) : length(first(points))
    if D == 2
        p = to_points(pts, 2)
        return PolygonDomain(convex_hull_2d(p))
    else
        p = to_points(pts, 3)
        normals, offsets = convex_hull_planes(p)
        return ConvexPolyhedronDomain(normals, offsets)
    end
end

function convex_hull_2d(p::Vector{Vec2})
    pts = sort(unique(p), by = x -> (x[1], x[2]))
    length(pts) < 3 && return pts
    lower = Vec2[]
    for x in pts
        while length(lower) >= 2 && cross2(lower[end] - lower[end - 1], x - lower[end]) <= 0
            pop!(lower)
        end
        push!(lower, x)
    end
    upper = Vec2[]
    for x in reverse(pts)
        while length(upper) >= 2 && cross2(upper[end] - upper[end - 1], x - upper[end]) <= 0
            pop!(upper)
        end
        push!(upper, x)
    end
    return vcat(lower[1:(end - 1)], upper[1:(end - 1)])
end

domain_dimension(::AbstractDomain{D}) where {D} = D

initial_polytope(d::PolygonDomain) = ClipPolygon(copy(d.vertices), [-k for k in 1:length(d.vertices)])
initial_polytope(d::ConvexPolyhedronDomain) = d.polytope

domain_bounding_box(d::PolygonDomain) = bounding_box(d.vertices)
domain_bounding_box(d::ConvexPolyhedronDomain) = bounding_box(d.polytope.verts)

function domain_diameter(d::AbstractDomain)
    lo, hi = domain_bounding_box(d)
    return norm(hi - lo)
end

"""
    inside(domain, x; tol = 0)

Test if `x` is inside the domain.
"""
inside(d::PolygonDomain, x; tol = 0.0) = inpolygon(Vec2(x), d.vertices; tol = tol)
function inside(d::ConvexPolyhedronDomain, x; tol = 0.0)
    for (n, o) in zip(d.normals, d.offsets)
        if dot(n, x) > o + tol
            return false
        end
    end
    return true
end

domain_volume(d::PolygonDomain) = abs(signed_area(d.vertices))
function domain_volume(d::ConvexPolyhedronDomain)
    P = d.polytope
    c = sum(P.verts) / length(P.verts)
    vol = 0.0
    for f in P.faces
        a = P.verts[f[1]]
        for i in 2:(length(f) - 1)
            vol += dot(cross(P.verts[f[i]] - a, P.verts[f[i + 1]] - a), a - c) / 6
        end
    end
    return vol
end
