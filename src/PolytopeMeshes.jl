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

end # module PolytopeMeshes
