"""
    PolytopeMeshes

A Julia package for conforming Voronoi meshes (Voronoi meshes) adapting to
geological features such as faults, fractures, and wells.

Based on the UPR module from MRST (MATLAB Reservoir Simulation Toolbox).

Reference:
  Berge, R.L., Klemetsdal, Ø.S. & Lie, K.-A. (2019).
  "Unstructured Voronoi meshes conforming to lower dimensional objects."
  Computational Geosciences, 23, 169–188.
"""
module PolytopeMeshes

using LinearAlgebra
using SparseArrays

# Mesh data structure
include("mesh.jl")

# Triangulation backend abstraction (before voronoi2D which uses it)
include("triangulation.jl")

# Utility functions
include("utils/geometry.jl")
include("utils/interpolation.jl")
include("utils/conflict.jl")
include("utils/mlqt.jl")

# 2D mesh generation
include("voronoi2D/clip_polygon.jl")
include("voronoi2D/split_at_intersections.jl")
include("voronoi2D/line_sites_2d.jl")
include("voronoi2D/surface_sites_2d.jl")
include("voronoi2D/surface_suf_cond_2d.jl")
include("voronoi2D/clipped_voronoi_2d.jl")
include("voronoi2D/sort_edges.jl")
include("voronoi2D/composite_voronoi_mesh_2d.jl")

# Exports
export UnstructuredMesh
export composite_voronoi_mesh_2d
export clipped_voronoi_2d
export line_sites_2d
export surface_sites_2d
export surface_suf_cond_2d
export remove_conflict_points
export compute_geometry!
export plot_mesh

# Triangulation backend exports
export TriangulationBackend
export BowyerWatson
export triangulate_2d
export default_triangulation_backend
export set_default_triangulation_backend!

"""
    plot_mesh(G::UnstructuredMesh; kwargs...)

Plot an `UnstructuredMesh`. Requires GLMakie to be loaded.

# Keyword Arguments
- `color_cells=true`: Color cells by index when `true`
- `colormap=:viridis`: Colormap for cell coloring
- `strokecolor=:black`: Color of cell edges
- `strokewidth=1.0`: Width of cell edges
- `show_nodes=false`: Show mesh nodes as scatter points
- `node_color=:red`: Color of node markers
- `node_size=4`: Size of node markers
- `figure_size=(800, 600)`: Size of the figure
- `title=""`: Title for the plot

Returns a `(fig, ax)` tuple.

!!! note
    This function requires the GLMakie package. Load it with `using GLMakie`
    before calling `plot_mesh`.
"""
function plot_mesh end

end # module PolytopeMeshes
