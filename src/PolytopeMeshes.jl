"""
    PolytopeMeshes

Generation of polytopal (Voronoi/PEBI) meshes in 2D and 3D that conform to
points, lines and surfaces. The site-placement algorithms are a port of the
UPR module of MRST (Berge, Klemetsdal & Lie, Comput Geosci 23, 2019).

Meshes are returned as [`PolytopalMesh`](@ref) and can be written to VTK with
[`write_vtk`](@ref), or converted to `Jutul.UnstructuredMesh` when Jutul is
loaded.
"""
module PolytopeMeshes

using LinearAlgebra
using Random
using StaticArrays
using NearestNeighbors
import DistMesh
import WriteVTK

const Vec2 = SVector{2, Float64}
const Vec3 = SVector{3, Float64}

include("mesh.jl")
include("geometry/points.jl")
include("geometry/polygons.jl")
include("geometry/polylines.jl")
include("geometry/planes.jl")
include("voronoi/clip2d.jl")
include("voronoi/clip3d.jl")
include("voronoi/domains.jl")
include("voronoi/voronoi.jl")
include("voronoi/centroidal.jl")
include("sites/conditions.jl")
include("sites/line_sites.jl")
include("sites/surface_sites_2d.jl")
include("sites/surface_sites_3d.jl")
include("triangulation/interface.jl")
include("triangulation/distmesh.jl")
include("sites/mlqt.jl")
include("grids/pebi_2d.jl")
include("grids/pebi_3d.jl")
include("grids/exact_3d.jl")
include("io/vtk.jl")
include("plotting.jl")

# Mesh type and queries
export PolytopalMesh
# number_of_cells/number_of_faces/number_of_nodes are not exported to avoid
# clashes with Jutul, which gets methods for PolytopalMesh in the extension.
export face_nodes, cell_faces, cell_nodes, is_boundary_face
export mesh_geometry, cell_volumes, cell_centroids, face_areas, face_normals, face_centroids
export check_mesh

# Domains and Voronoi diagrams
export PolygonDomain, ConvexPolyhedronDomain, box_domain, convex_hull_domain
export voronoi_mesh, centroidal_voronoi_mesh

# Site generation (UPR ports)
export split_at_intersections, interpolate_line_path
export line_sites_2d, line_sites_3d
export surface_sites_2d, surface_sites_3d, remove_surface_conflict_sites_3d
export surface_sufficient_condition, remove_conflict_points, line_sufficient_condition
export distmesh_2d, mlqt

# Triangulation backends
export TriangulationBackend, DistMeshBackend, triangulate_polygon

# Grid generators (UPR ports)
export clipped_pebi_grid_2d, pebi_grid_2d, composite_pebi_grid_2d
export composite_pebi_grid_3d, mirrored_pebi_grid_3d
export elliptic_surface_3d, surface_intersections_3d, line_grid_3d
export surface_grid_3d, volume_grid_3d, surface_triangulation

# Output
export write_vtk, write_vtk_faces

# Plotting (implemented in the Makie/GLMakie extensions)
export plot_polymesh, plot_polymesh!, plot_polymesh_data, plot_polymesh_data!, plot_polymesh_faces!, plot_polymesh_interactive

end # module PolytopeMeshes
