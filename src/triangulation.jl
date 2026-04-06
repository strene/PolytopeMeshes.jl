"""
    TriangulationBackend

Abstract type for 2D Delaunay triangulation backends.

Concrete subtypes must implement [`triangulate_2d`](@ref).
"""
abstract type TriangulationBackend end

"""
    BowyerWatson <: TriangulationBackend

Built-in Delaunay triangulation backend using the Bowyer-Watson algorithm.

This is the default backend when no external triangulation packages are loaded.
"""
struct BowyerWatson <: TriangulationBackend end

"""
    triangulate_2d(backend::TriangulationBackend, points::Matrix{Float64})

Compute 2D Delaunay triangulation of `points` using the given `backend`.

# Arguments
- `backend`: A [`TriangulationBackend`](@ref) instance
- `points`: n×2 matrix of point coordinates

# Returns
- `tri`: m×3 matrix of triangle vertex indices (each row is a triangle)
"""
function triangulate_2d end

# Module-level default backend
const _default_backend = Ref{TriangulationBackend}(BowyerWatson())

"""
    default_triangulation_backend()

Return the current default [`TriangulationBackend`](@ref).

When the `Gmsh` package is loaded, this automatically switches to
`GmshTriangulation`. Otherwise, it defaults to [`BowyerWatson`](@ref).
"""
function default_triangulation_backend()
    return _default_backend[]
end

"""
    set_default_triangulation_backend!(backend::TriangulationBackend)

Set the default [`TriangulationBackend`](@ref) used by mesh generation functions.
"""
function set_default_triangulation_backend!(backend::TriangulationBackend)
    _default_backend[] = backend
    return backend
end

"""
    generate_background_grid_2d(backend, fd, fh, h0, bbox, p_fix; kwargs...)

Generate a 2D background grid (point distribution and triangulation) using the
given `backend`.

When the backend is [`BowyerWatson`](@ref), this delegates to `distmesh_2d` for
iterative force-based point equilibration with Delaunay retriangulation.

Non-default backends (e.g. `GmshTriangulation`) may use their own mesh
generation strategies, bypassing DistMesh entirely.

# Arguments
- `backend`: A [`TriangulationBackend`](@ref) instance
- `fd`: Signed distance function `fd(p) -> Vector{Float64}`.
  Returns negative values inside the domain, positive outside.
- `fh`: Element size function `fh(p) -> Vector{Float64}`.
  Returns the desired relative edge length at each point.
- `h0`: Initial (base) edge length for the background grid.
- `bbox`: Bounding box `[x_min y_min; x_max y_max]` (2×2 matrix).
- `p_fix`: Fixed points that must appear in the output (n×2 matrix).
  These points appear first in the output.

# Keyword Arguments
- `max_iter`: Maximum number of iterations (default: 500, used by DistMesh)
- `poly_bdr`: Polygon boundary vertices (k×2 matrix, default: `zeros(0, 2)`).
  Used by backends that need explicit geometry (e.g. Gmsh).

# Returns
- `p`: Point coordinates (n×2 matrix). Fixed points are the first `nfix` rows.
- `t`: Triangulation (m×3 matrix of vertex indices).
"""
function generate_background_grid_2d end
