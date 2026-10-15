# Interchangeable triangulation backends. The grid generators only need a
# point set (and, for surface triangulations, triangles) in a polygon with a
# given size function and a set of points that must be included. Backends
# implement `triangulate_polygon` for their own subtype of
# `TriangulationBackend`, typically in a package extension.

"""
    TriangulationBackend

Abstract type for the backends used to distribute points and triangulate 2D
polygons, e.g. the background sites of [`voronoi_grid_2d`](@ref) and the
surface triangulations of [`surface_triangulation`](@ref). The default is
[`DistMeshBackend`](@ref).

A new backend is a subtype `B <: TriangulationBackend` with a method of
[`triangulate_polygon`](@ref). For a backend that depends on another package
(e.g. Gmsh.jl or Meshes.jl), define the type in PolytopeMeshes and the
method in a package extension, so that it becomes available when that
package is loaded.
"""
abstract type TriangulationBackend end

"""
    triangulate_polygon(backend, polygon, size, fixed; rng, min_size)

Triangulate the polygon with the given backend. Methods must accept:
- `polygon::Vector{SVector{2, Float64}}`: vertices of a simple polygon,
  ordered counter-clockwise.
- `size`: function `x -> h` giving the target edge length at `x`.
- `fixed::Vector{SVector{2, Float64}}`: distinct points in the closed polygon
  that must be vertices of the triangulation (constraint sites, polygon
  corners).
- `rng`: random number generator, for backends with random components.
- `min_size`: a lower bound of `size` in the polygon.

and return `(points, triangles)` where `points[1:length(fixed)]` are the
fixed points in the given order (up to round-off) and `triangles` is a
vector of `NTuple{3, Int}`. A backend that only distributes points may
return no triangles; it can then be used for the Voronoi grid generators but not
for [`surface_triangulation`](@ref). Points outside the polygon and degenerate
triangles are removed afterwards, so backends need not filter them.
"""
function triangulate_polygon(backend::TriangulationBackend, polygon, size, fixed; kwarg...)
    throw(ArgumentError(
        "triangulate_polygon is not implemented for $(typeof(backend)). " *
        "If the backend is provided by a package extension, load the package it depends on."
    ))
end

"""
    default_triangulation_backend()

The backend used when none is given: [`DistMeshBackend`](@ref)`()`.
"""
default_triangulation_backend() = DistMeshBackend()

# Call the backend, check that it honors the interface and clean up the result.
function _triangulate(backend::TriangulationBackend, polygon, size, fixed; rng = Random.default_rng(), min_size = nothing)
    polygon = to_points(polygon, 2)
    if signed_area(polygon) < 0
        polygon = reverse(polygon)
    end
    lo, hi = bounding_box(polygon)
    L = norm(hi - lo)
    fixed, _ = unique_points(to_points(fixed, 2); tol = 1e-10 * L)
    nfix = length(fixed)
    if isnothing(min_size)
        min_size = _min_size(x -> signed_distance_polygon(x, polygon), size, L / 100, lo, hi)
    end
    pts, tris = triangulate_polygon(backend, polygon, size, fixed; rng = rng, min_size = min_size)
    pts = Vec2.(pts)
    tris = [(Int(t[1]), Int(t[2]), Int(t[3])) for t in tris]
    if length(pts) < nfix || any(i -> norm(pts[i] - fixed[i]) > 1e-8 * L, 1:nfix)
        error("Backend $(typeof(backend)) did not return the fixed points first and in order")
    end
    pts[1:nfix] .= fixed
    # Remove free points outside the polygon, degenerate triangles and
    # triangles outside the polygon, then duplicate points.
    tol = 1e-3 * min_size
    keep = [i <= nfix || signed_distance_polygon(pts[i], polygon) < tol for i in eachindex(pts)]
    tris = filter(tris) do (i, j, k)
        keep[i] && keep[j] && keep[k] || return false
        area = abs(cross2(pts[j] - pts[i], pts[k] - pts[i])) / 2
        area > 1e-6 * min_size^2 && signed_distance_polygon((pts[i] + pts[j] + pts[k]) / 3, polygon) < 0
    end
    upts, map = unique_points(pts[keep]; tol = 1e-10 * L)
    newid = zeros(Int, length(pts))
    newid[keep] = map
    tris = [(newid[i], newid[j], newid[k]) for (i, j, k) in tris]
    tris = filter(t -> allunique(t), tris)
    return upts, tris, nfix
end
