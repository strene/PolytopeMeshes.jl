# Basic Voronoi Mesh Generation
#
# This example demonstrates how to create basic Voronoi meshes using
# composite_voronoi_mesh_2d, which is the Julia equivalent of MRST's
# compositePebiGrid2D.
#
# Based on: MRST book-ii/uprBookSection3.m (Section 3.1)

using PolytopeMeshes

## 1. Cartesian-equivalent Voronoi mesh
# The simplest use case: a rectangular domain with uniform cell sizes.
# This is equivalent to MRST's cartGrid, but produces a Voronoi mesh.
dx, dy = 0.25, 0.25   # mesh size in x and y direction
xmax, ymax = 1.0, 1.0 # domain dimensions

G, pts, F = composite_voronoi_mesh_2d([dx, dy], [xmax, ymax])
compute_geometry!(G)

println("=== Cartesian-equivalent Voronoi mesh ===")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Nodes: $(G.nodes.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println("  Domain: [0, $xmax] × [0, $ymax]")
println()

## 2. Finer mesh
# Reducing cell size produces more cells
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0])
compute_geometry!(G2)

println("=== Finer Voronoi mesh ===")
println("  Cells: $(G2.cells.num)")
println("  Faces: $(G2.faces.num)")
println("  Nodes: $(G2.nodes.num)")
println("  Total volume: $(round(sum(G2.cells.volumes), digits=4))")
println()

## 3. Rectangular domain with different aspect ratio
G3, _, _ = composite_voronoi_mesh_2d([0.15, 0.15], [2.0, 0.5])
compute_geometry!(G3)

println("=== Rectangular domain [2.0 × 0.5] ===")
println("  Cells: $(G3.cells.num)")
println("  Faces: $(G3.faces.num)")
println("  Total volume: $(round(sum(G3.cells.volumes), digits=4))")
