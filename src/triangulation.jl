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
