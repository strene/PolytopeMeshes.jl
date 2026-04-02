"""
    PolytopeMeshes

A Julia package for conforming Voronoi meshes (PEBI grids) adapting to
geological features such as faults, fractures, and wells.

Based on the UPR module from MRST (MATLAB Reservoir Simulation Toolbox).

Reference:
  Berge, R.L., Klemetsdal, Ø.S. & Lie, K.-A. (2019).
  "Unstructured Voronoi grids conforming to lower dimensional objects."
  Computational Geosciences, 23, 169–188.
"""
module PolytopeMeshes

using LinearAlgebra
using SparseArrays

# Grid data structure
include("grid.jl")

# Utility functions
include("utils/geometry.jl")
include("utils/interpolation.jl")
include("utils/conflict.jl")
include("utils/mlqt.jl")

# 2D grid generation
include("pebi2D/clip_polygon.jl")
include("pebi2D/split_at_intersections.jl")
include("pebi2D/line_sites_2d.jl")
include("pebi2D/surface_sites_2d.jl")
include("pebi2D/surface_suf_cond_2d.jl")
include("pebi2D/clipped_pebi_2d.jl")
include("pebi2D/sort_edges.jl")
include("pebi2D/composite_pebi_grid_2d.jl")

# Exports
export UnstructuredGrid
export composite_pebi_grid_2d
export clipped_pebi_2d
export line_sites_2d
export surface_sites_2d
export surface_suf_cond_2d
export remove_conflict_points
export compute_geometry!

end # module PolytopeMeshes
